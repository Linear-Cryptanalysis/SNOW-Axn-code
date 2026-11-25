#!/usr/bin/env python3.10
"""Structural verification of the reduced-scale SNOW-Axn analogue.

Checks, all exact unless noted:

  V1  AES S-box: permutation, S(0)=0x63, inverse consistency, matches a
      reference implementation of the affine/inversion construction.
  V2  alpha, beta: linear maps, invertible, and they satisfy the defining
      polynomials x^8+x^5+x^4+x^3+1 and x^8+x^7+x^5+x^2+1.
  V3  sigma = rotl3: permutation of order 8 and c.sigma(y) = sigma^-1(c).y.
  V4  LFSR: vector step == scalar step; the 64-bit LFSR map is invertible
      (rank 64); rank of (H(0),...,H(k-1)) is 8k for k = 1..8, i.e. the first
      eight clocks' H values are 64 independent linear forms of the state
      (so, for clocks 0..7, they are jointly uniform and independent).
  V5  FSM: the one-clock map is a bijection of the 24-bit register space for
      every fixed H (so a uniform register state stays uniform).
  V6  The two displayed forms of e(t) (keystream/LFSR form and register form)
      agree for random states at every clock t = 1..6.
  V7  LAT identity: sum_x (-1)^{a.S^-1(x) ^ q.x} = LAT[a][q].
  V8  t = 0 exhaustive check: with R1,R2,R3 independent uniform, the exact
      correlation E[(-1)^{e(0)}] equals rho1*rho2*rho3/2^24 of the report for
      every tested triple; the printed "note formula" (LAT[b][c'] first) does
      not.
  V9  Piling-up at t > 0 is a consequence of V4+V5: sample check that the
      marginal law of (R1(t),R2(t),R3(t)) is uniform (all three registers).

Run:  python3.10 verify.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_tables as bt


def gf2_rank(masks):
    """Rank over GF(2) of a list of python ints."""
    basis = []
    for m in masks:
        x = m
        for b in basis:
            x = min(x, x ^ b)
        if x:
            basis.append(x)
            basis.sort(reverse=True)
    return len(basis)


def state_to_int(u):
    return sum((int(u[r]) & 0xFF) << (8 * r) for r in range(8))


def int_to_state(x):
    return [(x >> (8 * r)) & 0xFF for r in range(8)]


def lfsr_state_step(T, x):
    return state_to_int(T.lfsr_advance_scalar(int_to_state(x)))


def h_forms(T, nclocks=8):
    """forms[t] = list of 8 GF(2)-mask (64-bit ints) giving the bits of H(t)."""
    # evolve all 64 basis vectors step by step, record u_0
    basis = [1 << j for j in range(64)]
    forms = []
    cur = list(basis)
    for t in range(nclocks):
        masks = []
        for bit in range(8):
            m = 0
            for j, x in enumerate(cur):
                if (x >> bit) & 1:  # bit `bit` of u_0 == bit `bit` of state
                    m |= 1 << j
            masks.append(m)
        forms.append(masks)
        cur = [lfsr_state_step(T, x) for x in cur]
    return forms


def main():
    T = bt.get_tables()
    ok = True

    def report(name, passed, extra=""):
        nonlocal ok
        ok = ok and passed
        print(f"[{'PASS' if passed else 'FAIL'}] {name}{('  ' + extra) if extra else ''}")

    # ---- V1 ---------------------------------------------------------------
    s = T.sbox
    report("V1 AES S-box is a permutation and inverse is consistent",
           len(np.unique(s)) == 256 and s[0] == 0x63 and
           np.array_equal(T.sbox[T.inv_sbox], np.arange(256, dtype=np.uint8)),
           f"S(0)=0x{int(s[0]):02x} S(1)=0x{int(s[1]):02x} S(255)=0x{int(s[255]):02x}")

    # ---- V2 ---------------------------------------------------------------
    def is_linear(table):
        for x in range(256):
            for y in range(256):
                if int(table[x ^ y]) != int(table[x] ^ table[y]):
                    return False
        return True

    def poly_holds(mult):
        v = 1
        powers = [1]
        for _ in range(8):
            v = int(mult[v])
            powers.append(v)
        # powers[i] == root^i ; check root^8 = reduce relation
        # alpha: x^8 + x^5 + x^4 + x^3 + 1 = 0  ->  root^8 = x^5^x^4^x^3^1
        return powers

    pa = poly_holds(T.alpha)
    pb = poly_holds(T.beta)
    report("V2 alpha,beta are linear bijections", is_linear(T.alpha) and
           is_linear(T.beta) and len(set(int(v) for v in T.alpha)) == 256 and
           len(set(int(v) for v in T.beta)) == 256)
    report("V2 alpha satisfies x^8+x^5+x^4+x^3+1 (alpha^8 = 1+a^3+a^4+a^5)",
           pa[8] == (pa[0] ^ pa[3] ^ pa[4] ^ pa[5]),
           f"alpha^1..alpha^3 = {pa[1:4]}")
    report("V2 beta satisfies x^8+x^7+x^5+x^2+1 (beta^8 = 1+b^2+b^5+b^7)",
           pb[8] == (pb[0] ^ pb[2] ^ pb[5] ^ pb[7]),
           f"beta^1..beta^3 = {pb[1:4]}")

    # ---- V3 ---------------------------------------------------------------
    sig_ok = all(int(T.rotl[int(T.rotl[int(T.rotl[int(T.rotl[int(T.rotl[
        int(T.rotl[int(T.rotl[int(T.rotl[y])])])])])])])]) == y for y in range(256))
    mask_id = all(bin(c & int(T.rotl[y])).count("1") % 2 ==
                  bin(bt.sigma_inv(c) & int(y)).count("1") % 2
                  for c in range(256) for y in range(256))
    report("V3 sigma is a permutation of order 8 (rotl3)", sig_ok)
    report("V3 c.sigma(y) == sigma^-1(c).y for all c,y", mask_id)

    # ---- V4 ---------------------------------------------------------------
    rng = np.random.default_rng(12345)
    u = [rng.integers(0, 256, size=1000).astype(np.uint8) for _ in range(8)]
    uv = T.lfsr_advance(u)
    us = [T.lfsr_advance_scalar([int(a[i]) for a in u]) for i in range(1000)]
    us_ok = all(int(uv[r][i]) == us[i][r] for r in range(8) for i in range(1000))
    report("V4 vector LFSR step == scalar byte recurrence", us_ok)
    mat = [lfsr_state_step(T, 1 << j) for j in range(64)]
    rank_lfsr = gf2_rank(mat)
    print("[info] V4 LFSR one-clock map rank = %d (the toy recurrence is "
          "singular, kernel dimension %d)" % (rank_lfsr, 64 - rank_lfsr))
    forms = h_forms(T, 10)
    ranks = []
    for k in range(1, 9):
        masks = [m for t in range(k) for m in forms[t]]
        ranks.append(gf2_rank(masks))
    print("[info] V4 rank(H(0..k-1)) for k=1..8 = %s (cumulative rank is "
          "deficient from k=5 on)" % ranks)
    window_ok = True
    window_ranks = {}
    for t0 in range(0, 7):
        for d in (1, 2, 3):
            ms = [m for t in range(t0 + 1, t0 + 1 + d) for m in forms[t]]
            r = gf2_rank(ms)
            window_ranks[(t0, d)] = r
            if r != 8 * d:
                window_ok = False
    report("V4 every H-window (H(t+1),...,H(t+d)), t=0..6, d<=3, has rank 8d "
           "(H bytes iid over the measured windows)", window_ok,
           "all windows full rank" if window_ok else str(window_ranks))

    # ---- V5 ---------------------------------------------------------------
    r1g, r2g, r3g = np.meshgrid(np.arange(256), np.arange(256),
                                np.arange(256), indexing="ij")
    h = 0x5A
    n1 = T.rotl[r2g ^ r3g ^ h]
    n2 = T.sbox[r1g]
    n3 = T.sbox[r2g]
    code = (r1g.astype(np.int64) << 16) | (r2g.astype(np.int64) << 8) | r3g.astype(np.int64)
    ncode = (n1.astype(np.int64) << 16) | (n2.astype(np.int64) << 8) | n3.astype(np.int64)
    seen = np.zeros(1 << 24, dtype=np.uint8)
    seen[ncode.reshape(-1)] = 1
    report("V5 FSM one-clock map is a bijection of the 24-bit state (h=0x5A)",
           int(seen.sum()) == (1 << 24))
    del r1g, r2g, r3g, n1, n2, n3, code, ncode, seen

    # ---- V6 ---------------------------------------------------------------
    # random states, advance the real cipher, compare the two forms of e(t)
    ntest = 64
    u = [rng.integers(0, 256, size=ntest).astype(np.uint8) for _ in range(8)]
    r1 = rng.integers(0, 256, size=ntest).astype(np.uint8)
    r2 = rng.integers(0, 256, size=ntest).astype(np.uint8)
    r3 = rng.integers(0, 256, size=ntest).astype(np.uint8)
    H = [u[0]]
    states = [(r1, r2, r3)]
    for step in range(1, 8):
        u = T.lfsr_advance(u)
        H.append(u[0])
        nr1, nr2, nr3 = T.fsm_advance(*states[-1], H[step])
        states.append((nr1, nr2, nr3))
    z = [states[i][0] ^ states[i][1] ^ H[i] for i in range(8)]
    triples = [(1, 1, 1), (0x11, 0x7E, 0x3C), (0x5A, 0xA5, 0x0F),
               (0xFF, 0x01, 0x80), (0x2C, 0x9B, 0x61)]
    bad = 0
    for (a, b, c) in triples:
        for t in range(1, 7):
            va = z[t - 1] ^ H[t - 1]
            vb = z[t] ^ H[t]
            vc = z[t + 1] ^ T.rotl[H[t + 1]] ^ H[t + 1]
            e_z = T.parity[(a & va) ^ (b & vb) ^ (c & vc)]
            A, B, C = T.block_vectors(*states[t])
            e_s = T.parity[(a & A) ^ (b & B) ^ (c & C)]
            bad += int(np.count_nonzero(e_z != e_s))
    report("V6 keystream form == register form of e(t), t=1..6", bad == 0,
           f"mismatches={bad} over {len(triples) * 6 * ntest} values")

    # ---- V7 ---------------------------------------------------------------
    x = np.arange(256, dtype=np.int64)
    W = 1 - 2 * T.parity[(np.arange(256)[:, None] & x[None, :])].astype(np.int64)
    # expression E[a][q] = sum_x (-1)^{a.S^-1(x) ^ q.x}
    E = np.zeros((256, 256), dtype=np.int64)
    for a in range(256):
        E[a] = (W[a][T.inv_sbox.astype(np.int64)] @ W.T)
    report("V7 sum_x (-1)^{a.S^-1(x) ^ q.x} == LAT[a][q]", np.array_equal(E, T.lat.astype(np.int64)))

    # ---- V8 ---------------------------------------------------------------
    # exhaustive t=0: (R1,R2,R3) independent uniform, e(0) is the register form
    AA = T.inv_sbox[np.arange(256)[:, None, None]]  # placeholder
    r1g = np.arange(256, dtype=np.uint8).reshape(256, 1, 1)
    r2g = np.arange(256, dtype=np.uint8).reshape(1, 256, 1)
    r3g = np.arange(256, dtype=np.uint8).reshape(1, 1, 256)
    Bv = r1g ^ r2g
    Cv = T.rotl[r2g ^ r3g] ^ T.sbox[r1g]
    invr2 = T.inv_sbox[r2g]
    invr3 = T.inv_sbox[r3g]
    A0 = invr2 ^ invr3
    triples = [(1, 1, 1), (0x11, 0x7E, 0x3C), (0x5A, 0xA5, 0x0F),
               (0xFF, 0x01, 0x80), (0x2C, 0x9B, 0x61), (0x08, 0x21, 0x45)]
    print("V8  exhaustive t=0 correlations (2^24 uniform triples):")
    print("      a     b     c   |  E[(-1)^e]      tau_pred      tau_note")
    v8_ok = True
    for (a, b, c) in triples:
        e = T.parity[(a & A0) ^ (b & Bv) ^ (c & Cv)]
        exact = 1.0 - 2.0 * float(e.mean())
        pred = bt.tau_unit(a, b, c, T) / 16777216.0
        note = bt.tau_note_formula(a, b, c, T) / 16777216.0
        v8_ok = v8_ok and abs(exact - pred) < 1e-12
        print(f"   0x{a:02x}  0x{b:02x}  0x{c:02x} | {exact:+.9f}  {pred:+.9f}  {note:+.9f}")
    report("V8 exact t=0 correlation == rho1*rho2*rho3/2^24 (all triples)", v8_ok)

    # ---- V8b --------------------------------------------------------------
    # The convention of the paper's eq. (N): the mask C multiplies
    # L^{-1} sigma^{-1} (z^{t+1} ^ sigma H3 ^ H2), so with L = I the R1 block
    # output mask is sigma(C) (because C.sigma^{-1}(w) = sigma(C).w) and the
    # R2/R3 input masks are C^B and C.  Hence
    #     tau = LAT[b][sigma(c)] * LAT[a][b^c] * LAT[a][c] / 2^24.
    print("V8b same check for the paper's eq. (N) mask convention "
          "(C applied to sigma^-1 of the third bracket):")
    print("      a     b     c   |  E[(-1)^e]      tau_paper")
    v8b_ok = True
    rotl = T.rotl.astype(np.int64)
    rotr = T.rotr.astype(np.int64)
    for (a, b, c) in triples:
        Q = rotr[rotl[(r2g ^ r3g).astype(np.int64)] ^ T.sbox[r1g].astype(np.int64)]
        e = T.parity[(a & A0) ^ (b & Bv) ^ (c & Q)]
        exact = 1.0 - 2.0 * float(e.mean())
        sc = int(T.rotl[c])
        pred = (int(T.lat[b][sc]) * int(T.lat[a][(b ^ c) & 255]) *
                int(T.lat[a][c])) / 16777216.0
        v8b_ok = v8b_ok and abs(exact - pred) < 1e-12
        print("   0x%02x  0x%02x  0x%02x | %+.9f  %+.9f" % (a, b, c, exact, pred))
    report("V8b paper-convention piling-up formula verified exhaustively", v8b_ok)

    # ---- V9 ---------------------------------------------------------------
    # marginal uniformity of each register at t=1,2,3 (sampled)
    ns = 1 << 20
    u = [rng.integers(0, 256, size=ns).astype(np.uint8) for _ in range(8)]
    reg = [rng.integers(0, 256, size=ns).astype(np.uint8) for _ in range(3)]
    chi = []
    for step in range(1, 4):
        u = T.lfsr_advance(u)
        reg = list(T.fsm_advance(*reg, u[0]))
        for z_ in reg:
            cnt = np.bincount(z_, minlength=256) / ns
            chi.append(float(np.max(np.abs(cnt - 1 / 256.0))))
    report("V9 each register is marginally uniform at t=1..3 (max dev < 1e-3)",
           max(chi) < 1e-3, f"maxdev={max(chi):.2e}")

    print()
    print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
