#!/usr/bin/env python3
"""Secure ECC TX -> AWGN channel -> RX BER/FER simulation.

Modes:
- uncoded: 640 information bits, no channel code.
- bch_only: BCH(840,640), hard-decision BCH decoder.
- ldpc_only: LDPC(1120,840), layered NMSA LDPC decoder.
- bch_ldpc: BCH(840,640) outer + LDPC(1120,840) inner.

This is a software/golden-model simulation for BER/FER plots. RTL simulation is
kept for functional verification and Vivado is used for resource/timing/power.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import random
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
PY_ROOT = THIS_DIR.parents[1]
LDPC_MODEL = PY_ROOT / "LDPC" / "ldpc_hex_model.py"


def import_ldpc_model():
    spec = importlib.util.spec_from_file_location("ldpc_hex_model", LDPC_MODEL)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


class BCHShortened:
    """Shortened primitive narrow-sense binary BCH over GF(2^10).

    The default is BCH(840,640,t=20), shortened from length 1023. It is meant as
    the software golden model matching the current project choice of 5 AES
    blocks = 640 information bits -> 840 BCH codeword bits.
    """

    def __init__(self, n: int = 840, k: int = 640, t: int = 20) -> None:
        self.m = 10
        self.n_full = (1 << self.m) - 1
        self.n = n
        self.k = k
        self.r = n - k
        self.t = t
        self.prim_poly = 0x409
        self.gf_exp, self.gf_log = self._make_gf()
        self.g = self._generator_poly()
        self.g_degree = len(self.g) - 1
        if self.g_degree > self.r:
            raise RuntimeError(
                f"BCH generator degree={self.g_degree}, expected parity={self.r}. "
                "Choose a BCH(n,k,t) combination with enough parity bits."
            )
        self.pad_parity = self.r - self.g_degree

    def _make_gf(self) -> tuple[list[int], list[int]]:
        gf_exp = [0] * (2 * self.n_full)
        gf_log = [-1] * (self.n_full + 1)
        x = 1
        for i in range(self.n_full):
            gf_exp[i] = x
            gf_log[x] = i
            x <<= 1
            if x & (1 << self.m):
                x ^= self.prim_poly
        for i in range(self.n_full, 2 * self.n_full):
            gf_exp[i] = gf_exp[i - self.n_full]
        return gf_exp, gf_log

    def _gf_mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.gf_exp[self.gf_log[a] + self.gf_log[b]]

    def _gf_div(self, a: int, b: int) -> int:
        if a == 0:
            return 0
        if b == 0:
            raise ZeroDivisionError("GF division by zero")
        return self.gf_exp[(self.gf_log[a] - self.gf_log[b]) % self.n_full]

    def _minimal_poly(self, i: int) -> list[int]:
        seen: list[int] = []
        x = i % self.n_full
        while x not in seen:
            seen.append(x)
            x = (2 * x) % self.n_full
        poly = [1]
        for pwr in seen:
            root = self.gf_exp[pwr]
            new = [0] * (len(poly) + 1)
            for j, coef in enumerate(poly):
                new[j] ^= self._gf_mul(coef, root)
                new[j + 1] ^= coef
            poly = new
        return [c & 1 for c in poly]

    @staticmethod
    def _poly_mul_bin(a: list[int], b: list[int]) -> list[int]:
        out = [0] * (len(a) + len(b) - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        out[i + j] ^= 1
        return out

    def _generator_poly(self) -> list[int]:
        used = set()
        g = [1]
        for i in range(1, 2 * self.t + 1):
            cyc: list[int] = []
            x = i % self.n_full
            while x not in cyc:
                cyc.append(x)
                x = (2 * x) % self.n_full
            key = min(cyc)
            if key in used:
                continue
            used.add(key)
            g = self._poly_mul_bin(g, self._minimal_poly(i))
        return g

    def encode_bits(self, msg: list[int]) -> list[int]:
        if len(msg) != self.k:
            raise ValueError(f"BCH message must have {self.k} bits")
        work = [0] * (self.k + self.g_degree)
        work[self.g_degree:] = [b & 1 for b in msg]
        for i in range(self.k - 1, -1, -1):
            if work[self.g_degree + i]:
                for j, gj in enumerate(self.g):
                    if gj:
                        work[i + j] ^= 1
        parity = work[: self.g_degree] + [0] * self.pad_parity
        return [b & 1 for b in msg] + parity

    def _syndromes_full(self, rx: list[int]) -> list[int]:
        full = [0] * self.n_full
        active_parity = rx[self.k : self.k + self.g_degree]
        full[: self.g_degree] = active_parity
        full[self.g_degree : self.g_degree + self.k] = rx[: self.k]
        ones = [i for i, bit in enumerate(full) if bit]
        syndromes: list[int] = []
        for j in range(1, 2 * self.t + 1):
            s = 0
            for pos in ones:
                s ^= self.gf_exp[(j * pos) % self.n_full]
            syndromes.append(s)
        return syndromes

    def _berlekamp_massey(self, synd: list[int]) -> tuple[list[int], int]:
        c = [1] + [0] * (2 * self.t)
        b = [1] + [0] * (2 * self.t)
        l_val = 0
        m_val = 1
        bb = 1
        for n_idx in range(2 * self.t):
            d = synd[n_idx]
            for i in range(1, l_val + 1):
                d ^= self._gf_mul(c[i], synd[n_idx - i])
            if d == 0:
                m_val += 1
                continue
            old = c[:]
            coef = self._gf_div(d, bb)
            for i in range(0, 2 * self.t + 1 - m_val):
                if b[i]:
                    c[i + m_val] ^= self._gf_mul(coef, b[i])
            if 2 * l_val <= n_idx:
                l_val = n_idx + 1 - l_val
                b = old
                bb = d
                m_val = 1
            else:
                m_val += 1
        return c[: l_val + 1], l_val

    def decode_codeword(self, rx: list[int]) -> tuple[list[int], list[int], bool, int]:
        """Decode a shortened BCH word and retain its corrected codeword.

        The returned codeword uses the same systematic order as ``encode_bits``:
        message bits followed by parity bits. Keeping this representation is
        needed by BCH-aided LDPC studies, where all corrected BCH bits,
        including parity bits, are fed back to the inner decoder.
        """
        if len(rx) != self.n:
            raise ValueError(f"BCH codeword must have {self.n} bits")
        received = [bit & 1 for bit in rx]
        full = [0] * self.n_full
        active_parity = received[self.k : self.k + self.g_degree]
        full[: self.g_degree] = active_parity
        full[self.g_degree : self.g_degree + self.k] = received[: self.k]
        synd = self._syndromes_full(received)
        if max(synd) == 0:
            return received[: self.k], received, True, 0
        sigma, degree = self._berlekamp_massey(synd)
        error_positions: list[int] = []
        for pos in range(self.n_full):
            x_inv = self.gf_exp[(self.n_full - pos) % self.n_full]
            val = 0
            xp = 1
            for coef in sigma:
                if coef:
                    val ^= self._gf_mul(coef, xp)
                xp = self._gf_mul(xp, x_inv)
            if val == 0:
                error_positions.append(pos)
        ok = len(error_positions) == degree
        for pos in error_positions:
            if pos >= self.g_degree + self.k:
                ok = False
            else:
                full[pos] ^= 1
        corrected = full[self.g_degree : self.g_degree + self.k]
        check_rx = corrected + full[: self.g_degree] + [0] * self.pad_parity
        ok = ok and max(self._syndromes_full(check_rx)) == 0
        return corrected, check_rx, ok, len(error_positions)

    def decode_bits(self, rx: list[int]) -> tuple[list[int], bool, int]:
        """Backward-compatible payload-only BCH decode API."""
        corrected, _codeword, ok, correction_count = self.decode_codeword(rx)
        return corrected, ok, correction_count


def rand_bits(rng: random.Random, n_bits: int) -> list[int]:
    return [rng.getrandbits(1) for _ in range(n_bits)]


def count_errors(a: list[int], b: list[int]) -> int:
    return sum((x ^ y) & 1 for x, y in zip(a, b))


def awgn_bpsk_llr(bits: list[int], ebn0_db: float, rate: float, rng: random.Random) -> list[float]:
    ebn0 = 10.0 ** (ebn0_db / 10.0)
    sigma = math.sqrt(1.0 / (2.0 * rate * ebn0))
    inv_var2 = 2.0 / (sigma * sigma)
    llr: list[float] = []
    for bit in bits:
        x = 1.0 if bit == 0 else -1.0
        y = x + rng.gauss(0.0, sigma)
        llr.append(inv_var2 * y)
    return llr


def hard_from_llr(llr: list[float]) -> list[int]:
    return [0 if value >= 0.0 else 1 for value in llr]


def simulate_frame(mode: str, ebn0_db: float, rng: random.Random, ldpc, bch: BCHShortened, llr_mag_clip: float | None) -> tuple[int, int, bool, int, int]:
    if mode == "uncoded":
        info = rand_bits(rng, 640)
        llr = awgn_bpsk_llr(info, ebn0_db, 1.0, rng)
        dec = hard_from_llr(llr)
        bit_err = count_errors(info, dec)
        return bit_err, len(info), bit_err > 0, 0, 0

    if mode == "bch_only":
        info = rand_bits(rng, 640)
        code = bch.encode_bits(info)
        llr = awgn_bpsk_llr(code, ebn0_db, 640 / 840, rng)
        hard = hard_from_llr(llr)
        dec, ok, corrected = bch.decode_bits(hard)
        bit_err = count_errors(info, dec)
        return bit_err, len(info), bit_err > 0 or not ok, 0, 0 if ok else 1

    if mode == "ldpc_only":
        info = rand_bits(rng, 840)
        code = ldpc.encode_bits(info)
        llr = awgn_bpsk_llr(code, ebn0_db, 840 / 1120, rng)
        if llr_mag_clip is not None:
            llr = [max(-llr_mag_clip, min(llr_mag_clip, x)) for x in llr]
        dec, hard, ok, it_used = ldpc.decode_llr(llr)
        bit_err = count_errors(info, dec)
        return bit_err, len(info), bit_err > 0 or not ok, it_used, 0 if ok else 1

    if mode == "bch_ldpc":
        info = rand_bits(rng, 640)
        bch_code = bch.encode_bits(info)
        ldpc_code = ldpc.encode_bits(bch_code)
        llr = awgn_bpsk_llr(ldpc_code, ebn0_db, 640 / 1120, rng)
        if llr_mag_clip is not None:
            llr = [max(-llr_mag_clip, min(llr_mag_clip, x)) for x in llr]
        ldpc_dec, hard, ldpc_ok, it_used = ldpc.decode_llr(llr)
        bch_dec, bch_ok, corrected = bch.decode_bits(ldpc_dec)
        bit_err = count_errors(info, bch_dec)
        frame_err = bit_err > 0 or not ldpc_ok or not bch_ok
        return bit_err, len(info), frame_err, it_used, (0 if ldpc_ok else 1) + (0 if bch_ok else 1)

    raise ValueError(f"Unsupported mode {mode}")


def parse_snr_list(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(description="Secure ECC BER/FER link simulation")
    parser.add_argument("--snr", default="0,1,2,3,4,5", help="Comma-separated Eb/N0 values in dB")
    parser.add_argument("--frames", type=int, default=100, help="Frames per SNR per mode")
    parser.add_argument("--max-frame-errors", type=int, default=50, help="Stop a point early after this many frame errors; 0 disables")
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--max-iter", type=int, default=30)
    parser.add_argument("--alpha", type=float, default=0.75)
    parser.add_argument("--llr-clip", type=float, default=31.0, help="Clip LDPC LLRs to mimic fixed-point range; <=0 disables")
    parser.add_argument("--modes", default="uncoded,bch_only,ldpc_only,bch_ldpc")
    parser.add_argument("--out", type=Path, default=THIS_DIR / "results" / "ber_fer_results.csv")
    args = parser.parse_args()

    ldpc_mod = import_ldpc_model()
    cfg = dict(ldpc_mod.PRESETS["1120_840"])
    cfg["alpha"] = args.alpha
    cfg["max_iter"] = args.max_iter
    ldpc = ldpc_mod.QCLDPC(cfg["pcm"], int(cfg["z"]), float(cfg["alpha"]), int(cfg["max_iter"]))
    bch = BCHShortened(n=840, k=640, t=20)

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    snrs = parse_snr_list(args.snr)
    rng = random.Random(args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    llr_clip = None if args.llr_clip <= 0 else args.llr_clip

    rows: list[dict[str, object]] = []
    print(f"LDPC: {ldpc.summary()}")
    print(f"BCH : N={bch.n}, K={bch.k}, t={bch.t}, parity={bch.r}, active_parity={bch.g_degree}, reserved={bch.pad_parity}")
    print(f"Output CSV: {args.out}")

    for mode in modes:
        for snr in snrs:
            bit_errors = 0
            total_bits = 0
            frame_errors = 0
            total_iter = 0
            decode_fail = 0
            frames_done = 0
            for _ in range(args.frames):
                be, nb, fe, it, df = simulate_frame(mode, snr, rng, ldpc, bch, llr_clip)
                bit_errors += be
                total_bits += nb
                frame_errors += int(fe)
                total_iter += it
                decode_fail += df
                frames_done += 1
                if args.max_frame_errors > 0 and frame_errors >= args.max_frame_errors:
                    break
            ber = bit_errors / total_bits if total_bits else 0.0
            fer = frame_errors / frames_done if frames_done else 0.0
            avg_iter = total_iter / frames_done if frames_done else 0.0
            row = {
                "mode": mode,
                "ebn0_db": snr,
                "frames": frames_done,
                "bits": total_bits,
                "bit_errors": bit_errors,
                "frame_errors": frame_errors,
                "ber": ber,
                "fer": fer,
                "avg_iter": avg_iter,
                "decode_fail": decode_fail,
            }
            rows.append(row)
            print(
                f"{mode:9s} Eb/N0={snr:4.1f} dB frames={frames_done:5d} "
                f"BER={ber:.3e} FER={fer:.3e} avg_iter={avg_iter:.2f} fail={decode_fail}"
            )

    with args.out.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
