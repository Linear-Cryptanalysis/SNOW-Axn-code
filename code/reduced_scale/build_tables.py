#!/usr/bin/env python3.10
"""Tables for the reduced-scale SNOW-Axn analogue.

The toy cipher (all parameters fixed by ``_rev/fse2`` work order):

  * word size 8 bit, FSM registers R1,R2,R3 in GF(2)^8,
  * S = AES S-box (computed from scratch, no external data),
  * sigma = cyclic left rotation of the byte by 3 bits (a linear permutation),
  * E = S because the linear layer is L = identity,
  * LFSR: eight bytes u_0..u_7 with the byte recurrence of the design
        u'_{r} = c_r*u_{r} ^ u_{r+1} ^ u_{r+3},  r = 0..7  (indices mod 8),
    where c_r = alpha for r < 4 and c_r = beta for r >= 4.  alpha is a root of
    x^8 + x^5 + x^4 + x^3 + 1 and beta a root of x^8 + x^7 + x^5 + x^2 + 1, so
    multiplication by alpha / beta is the companion-matrix map of that
    polynomial in the polynomial basis 1, alpha, ..., alpha^7.  Only those two
    linear maps are ever used, so the toy does not need a full GF(2^8)
    multiplication table for alpha/beta (the AES field is used for the S-box).

  * H(t) = u_0(t),  z(t) = R1(t) ^ R2(t) ^ H(t),
    R1(t+1) = sigma(R2(t) ^ R3(t) ^ H(t+1)),
    R2(t+1) = S(R1(t)),  R3(t+1) = S(R2(t)).

Linear approximation (three clocks, as in the paper's eq. (N) with L = I):

  e(t) = a.[z(t-1) ^ H(t-1)] ^ b.[z(t) ^ H(t)]
         ^ c.[z(t+1) ^ sigma H(t+1) ^ H(t+1)]
       = a.[S^-1(R2) ^ S^-1(R3)] ^ b.[R1 ^ R2] ^ c.[sigma(R2 ^ R3) ^ S(R1)]

(second line at clock t; verify.py checks the equality numerically).  Because
sigma is a bit permutation,  c . sigma(y) = sigma^-1(c) . y ; writing
c' = sigma^-1(c) the approximation splits into three one-register blocks

  e = e_R1(R1) ^ e_R2(R2) ^ e_R3(R3),
  e_R1(x) = b.x ^ c.S(x),
  e_R2(x) = a.S^-1(x) ^ (b ^ c').x,
  e_R3(x) = a.S^-1(x) ^ c'.x,

so that with R1,R2,R3 independent and uniform

  tau = rho1 * rho2 * rho3 / 2^24,
  rho1 = LAT[b][c],  rho2 = LAT[a][b^c'],  rho3 = LAT[a][c'],

where LAT[a][b] = sum_x (-1)^{a.x ^ b.S(x)}.  (The identity
sum_x (-1)^{a.S^-1(x) ^ q.x} = LAT[a][q] is used for the two inverse blocks;
verify.py checks it and the exhaustive t = 0 correlation.)  Note that this is
the piling-up value of the toy *as defined above*: the first factor is
LAT[b][c], not LAT[b][c'] -- see verify.py ``--note-formula`` and the report.

All tables are returned as numpy arrays; the module has no side effects other
than building them on import (cheap, milliseconds).
"""

import numpy as np

MASK = 0xFF
SHIFT = 3  # sigma = rotl3

# ---------------------------------------------------------------------------
# GF(2^8) for the AES S-box (polynomial x^8+x^4+x^3+x+1 = 0x11B)
# ---------------------------------------------------------------------------

def gf_mul(a, b, poly=0x11B):
    """Product in GF(2^8) modulo ``poly`` (slow reference implementation)."""
    r = 0
    for _ in range(8):
        if b & 1:
            r ^= a
        b >>= 1
        hi = a & 0x80
        a = (a << 1) & MASK
        if hi:
            a ^= poly & MASK
    return r


def _parity_table():
    return np.array([bin(x).count("1") & 1 for x in range(256)], dtype=np.uint8)


def build_aes_sbox():
    """AES S-box, computed from scratch: inversion in GF(2^8) + affine map."""
    inv = np.zeros(256, dtype=np.uint8)
    for x in range(1, 256):
        for y in range(1, 256):
            if gf_mul(x, y) == 1:
                inv[x] = y
                break
    def rotl8(v, s):
        return ((v << s) | (v >> (8 - s))) & MASK
    sbox = np.zeros(256, dtype=np.uint8)
    for x in range(256):
        v = int(inv[x])
        sbox[x] = v ^ rotl8(v, 1) ^ rotl8(v, 2) ^ rotl8(v, 3) ^ rotl8(v, 4) ^ 0x63
    return sbox


# ---------------------------------------------------------------------------
# alpha / beta : companion-matrix linear maps of the two design polynomials
# ---------------------------------------------------------------------------

ALPHA_REDUCE = 0x39  # alpha^8 = 1 + alpha^3 + alpha^4 + alpha^5
BETA_REDUCE = 0xA5   # beta^8  = 1 + beta^2 + beta^5 + beta^7


def _mul_by_root_table(reduce_byte):
    """Linear map table v -> root * v in the polynomial basis 1, root, ..."""
    t = np.zeros(256, dtype=np.uint8)
    for v in range(256):
        out = (v << 1) & MASK
        if v & 0x80:
            out ^= reduce_byte
        t[v] = out
    return t


def _rotl_table(s):
    return np.array([((v << s) | (v >> (8 - s))) & MASK for v in range(256)],
                    dtype=np.uint8)


# ---------------------------------------------------------------------------
# LAT
# ---------------------------------------------------------------------------

def build_lat(sbox):
    """LAT[a][b] = sum_x (-1)^{a.x ^ b.S(x)} as int16 (range [-256,256])."""
    x = np.arange(256, dtype=np.int64)
    par = _parity_table()
    # W[a][x] = (-1)^{a.x}
    W = 1 - 2 * par[(np.arange(256)[:, None] & x[None, :])].astype(np.int16)
    # Ws[b][x] = (-1)^{b.S(x)}
    Ws = W[:, sbox.astype(np.int64)]
    return (W.astype(np.int32) @ Ws.T.astype(np.int32)).astype(np.int16)


# ---------------------------------------------------------------------------
# assembled table bundle
# ---------------------------------------------------------------------------

class Tables:
    __slots__ = ("sbox", "inv_sbox", "parity", "lat", "alpha", "beta",
                 "rotl", "rotr")

    def __init__(self):
        self.sbox = build_aes_sbox()
        self.inv_sbox = np.zeros(256, dtype=np.uint8)
        self.inv_sbox[self.sbox] = np.arange(256, dtype=np.uint8)
        self.parity = _parity_table()
        self.lat = build_lat(self.sbox)
        self.alpha = _mul_by_root_table(ALPHA_REDUCE)
        self.beta = _mul_by_root_table(BETA_REDUCE)
        self.rotl = _rotl_table(SHIFT)
        self.rotr = _rotl_table(8 - SHIFT)

    # -- scalar cipher primitives (used by verify.py / exact_folded.py) ----
    def lfsr_advance_scalar(self, u):
        a, b = int(self.alpha[0]), int(self.beta[0])  # unused, keep explicit
        u = [int(v) & MASK for v in u]
        return [
            int(self.alpha[u[0]]) ^ u[1] ^ u[3],
            int(self.alpha[u[1]]) ^ u[2] ^ u[4],
            int(self.alpha[u[2]]) ^ u[3] ^ u[5],
            int(self.alpha[u[3]]) ^ u[4] ^ u[6],
            int(self.beta[u[4]]) ^ u[5] ^ u[7],
            int(self.beta[u[5]]) ^ u[6] ^ u[0],
            int(self.beta[u[6]]) ^ u[7] ^ u[1],
            int(self.beta[u[7]]) ^ u[0] ^ u[2],
        ]

    # -- vectorised LFSR step, u = list of 8 uint8 arrays -------------------
    def lfsr_advance(self, u):
        u0, u1, u2, u3, u4, u5, u6, u7 = u
        a, b = self.alpha, self.beta
        return [a[u0] ^ u1 ^ u3,
                a[u1] ^ u2 ^ u4,
                a[u2] ^ u3 ^ u5,
                a[u3] ^ u4 ^ u6,
                b[u4] ^ u5 ^ u7,
                b[u5] ^ u6 ^ u0,
                b[u6] ^ u7 ^ u1,
                b[u7] ^ u0 ^ u2]

    # -- vectorised FSM step, registers uint8 arrays ------------------------
    def fsm_advance(self, r1, r2, r3, h_next):
        n1 = self.rotl[r2 ^ r3 ^ h_next]
        n2 = self.sbox[r1]
        n3 = self.sbox[r2]
        return n1, n2, n3

    # -- linear-approximation block vectors (uint8 arrays) ------------------
    def block_vectors(self, r1, r2, r3):
        """Return the three byte-valued block quantities A,B,C with
        e = a.A ^ b.B ^ c.C for any linear masks a,b,c."""
        A = self.inv_sbox[r2] ^ self.inv_sbox[r3]
        B = r1 ^ r2
        C = self.rotl[r2 ^ r3] ^ self.sbox[r1]
        return A, B, C


TABLES = None


def get_tables():
    global TABLES
    if TABLES is None:
        TABLES = Tables()
    return TABLES


# ---------------------------------------------------------------------------
# prediction helpers
# ---------------------------------------------------------------------------

def sigma_inv(c):
    """c' = sigma^-1(c) (rotr by 3) for a 0..255 integer."""
    return int(((c >> SHIFT) | (c << (8 - SHIFT))) & MASK)


def tau_unit(a, b, c, T=None):
    """Integer numerator of the piling-up correlation: rho1*rho2*rho3.
    tau = tau_unit / 2^24."""
    T = T or get_tables()
    cp = sigma_inv(c)
    return int(T.lat[b][c]) * int(T.lat[a][(b ^ cp) & MASK]) * int(T.lat[a][cp])


def tau(a, b, c, T=None):
    return tau_unit(a, b, c, T) / 16777216.0


def tau_note_formula(a, b, c, T=None):
    """The expression as printed in the work order, i.e. with LAT[b][c'] as the
    first factor.  Only used by verify.py to document that it does NOT match
    the toy as defined; kept for the report."""
    T = T or get_tables()
    cp = sigma_inv(c)
    return int(T.lat[b][cp]) * int(T.lat[a][(b ^ cp) & MASK]) * int(T.lat[a][cp])


if __name__ == "__main__":
    T = get_tables()
    print("AES S-box [0:16] =", [int(v) for v in T.sbox[:16]])
    print("alpha table [0:8] =", [int(v) for v in T.alpha[:8]])
    print("beta  table [0:8] =", [int(v) for v in T.beta[:8]])
    print("sigma identity: c.sigma(y) == sigma^-1(c).y for all c,y:",
          all(bin(c & int(T.rotl[y])).count("1") % 2 ==
              bin(sigma_inv(c) & int(y)).count("1") % 2
              for c in range(256) for y in range(256)))
    print("max |LAT| =", int(np.abs(T.lat).max()))
