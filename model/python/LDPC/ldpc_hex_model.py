#!/usr/bin/env python3
"""Generic QC-LDPC hex-input/hex-output model.

Features:
- Configurable Z, PCM/base matrix, alpha, max_iter.
- Built-in presets:
  - 1120_840: Z=56, PCM 5x20
  - 128_64  : Z=16, PCM 4x8
- Systematic encoding using H = [A | B], p = B^-1 A m over GF(2).
- Layered normalized min-sum decoding from hard codeword hex mapped to LLR.

Default padded mode:
- encode input : payload hex bytes
- internal     : zero-pad to K-bit blocks, no length header
- encode output: LDPC codeword hex
- decode output: full recovered payload hex, including padding if any

Raw block mode:
- use --raw-blocks
- encode input length must be a multiple of K bits and byte-aligned
- decode input length must be a multiple of N bits and byte-aligned
"""

from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

PRESETS = {
    "1120_840": {
        "z": 56,
        "alpha": 0.75,
        "max_iter": 50,
        "pcm": [
            [5, 14, 12, 1, 2, 37, 45, 26, 24, 0, 3, -1, 34, 7, 46, 10, -1, -1, -1, -1],
            [0, 35, 1, 26, 0, 10, 16, 16, 34, 4, 2, 23, 0, 51, -1, 49, 20, -1, -1, -1],
            [12, 28, 22, 46, 3, 16, 51, 2, 25, 29, 19, 18, 52, -1, 37, -1, 34, 39, -1, -1],
            [0, 51, 16, 31, 13, 39, 27, 33, 8, 27, 53, 13, -1, 52, 33, -1, -1, 38, 7, -1],
            [36, 6, 3, 51, 4, 19, 4, 45, 48, 9, -1, 11, 22, 23, 43, -1, -1, -1, 14, 1],
        ],
    },
    "128_64": {
        "z": 16,
        "alpha": 0.75,
        "max_iter": 30,
        "pcm": [
            [3, 8, 14, -1, 0, -1, -1, -1],
            [9, -1, 1, 5, 0, 0, -1, -1],
            [2, 11, -1, 7, -1, 0, 0, -1],
            [-1, 4, 10, 13, 0, -1, 0, 0],
        ],
    },
}


def read_hex_file(path: Path) -> bytes:
    text = "".join(path.read_text(encoding="utf-8-sig").split())
    if text.startswith("0x") or text.startswith("0X"):
        text = text[2:]
    if len(text) % 2:
        raise ValueError("Hex input must contain an even number of hex digits")
    return bytes.fromhex(text)


def write_hex_file(path: Path, data: bytes) -> None:
    path.write_text(data.hex(), encoding="utf-8")


def bytes_to_bits(data: bytes) -> list[int]:
    return [(byte >> bit) & 1 for byte in data for bit in range(7, -1, -1)]


def bits_to_bytes(bits: list[int]) -> bytes:
    if len(bits) % 8 != 0:
        raise ValueError("Number of bits must be divisible by 8")
    out = bytearray()
    for i in range(0, len(bits), 8):
        value = 0
        for bit in bits[i : i + 8]:
            value = (value << 1) | (bit & 1)
        out.append(value)
    return bytes(out)


def parse_positions(text: str | None) -> list[int]:
    if not text:
        return []
    out: list[int] = []
    for item in text.replace(";", ",").split(","):
        item = item.strip()
        if item:
            out.append(int(item, 0))
    return out


def load_config(args: argparse.Namespace) -> dict:
    if args.config:
        cfg = json.loads(args.config.read_text(encoding="utf-8-sig"))
    else:
        cfg = dict(PRESETS[args.preset])
    if args.z is not None:
        cfg["z"] = args.z
    if args.alpha is not None:
        cfg["alpha"] = args.alpha
    if args.max_iter is not None:
        cfg["max_iter"] = args.max_iter
    return cfg


class QCLDPC:
    def __init__(self, pcm: list[list[int]], z: int, alpha: float = 0.75, max_iter: int = 50) -> None:
        self.pcm = pcm
        self.z = z
        self.alpha = alpha
        self.max_iter = max_iter
        self.mb = len(pcm)
        self.nb = len(pcm[0])
        if any(len(row) != self.nb for row in pcm):
            raise ValueError("PCM rows must have the same length")
        self.n = self.nb * z
        self.m = self.mb * z
        self.k = self.n - self.m
        if self.k <= 0:
            raise ValueError("Invalid PCM/Z: K=N-M must be positive")
        if self.n % 8 or self.k % 8:
            raise ValueError("This file model requires byte-aligned N and K")
        self.info_bytes = self.k // 8
        self.codeword_bytes = self.n // 8
        self.layers = [[(j, row[j]) for j in range(self.nb) if row[j] >= 0] for row in pcm]
        self.h_rows = self._build_h_rows()
        self.a_rows = [row & ((1 << self.k) - 1) for row in self.h_rows]
        self.b_rows = [row >> self.k for row in self.h_rows]
        self.b_inv_rows = self._gf2_inverse_rows(self.b_rows, self.m)

    def _build_h_rows(self) -> list[int]:
        rows = [0] * self.m
        for bi in range(self.mb):
            for bj in range(self.nb):
                sh = self.pcm[bi][bj]
                if sh < 0:
                    continue
                for r in range(self.z):
                    col = bj * self.z + ((r + sh) % self.z)
                    rows[bi * self.z + r] |= 1 << col
        return rows

    @staticmethod
    def _gf2_inverse_rows(rows: list[int], n: int) -> list[int]:
        aug = [rows[i] | (1 << (n + i)) for i in range(n)]
        pivot_row = 0
        for col in range(n):
            pivot = None
            for r in range(pivot_row, n):
                if (aug[r] >> col) & 1:
                    pivot = r
                    break
            if pivot is None:
                continue
            aug[pivot_row], aug[pivot] = aug[pivot], aug[pivot_row]
            for r in range(n):
                if r != pivot_row and ((aug[r] >> col) & 1):
                    aug[r] ^= aug[pivot_row]
            pivot_row += 1
        if pivot_row != n:
            raise RuntimeError("Parity submatrix B is singular. The last M bits cannot be used directly as parity.")
        mask = (1 << n) - 1
        return [(row >> n) & mask for row in aug]

    @staticmethod
    def _bits_to_int(bits: list[int]) -> int:
        out = 0
        for i, bit in enumerate(bits):
            if bit & 1:
                out |= 1 << i
        return out

    @staticmethod
    def _int_to_bits(value: int, width: int) -> list[int]:
        return [(value >> i) & 1 for i in range(width)]

    @staticmethod
    def _parity(x: int) -> int:
        return x.bit_count() & 1

    def encode_bits(self, msg_bits: list[int]) -> list[int]:
        if len(msg_bits) != self.k:
            raise ValueError(f"LDPC message must have {self.k} bits")
        msg_int = self._bits_to_int(msg_bits)
        syndrome_bits = [self._parity(row & msg_int) for row in self.a_rows]
        syndrome_int = self._bits_to_int(syndrome_bits)
        parity_bits = [self._parity(inv_row & syndrome_int) for inv_row in self.b_inv_rows]
        return msg_bits + parity_bits

    def check_bits(self, bits: list[int]) -> bool:
        if len(bits) != self.n:
            return False
        value = self._bits_to_int(bits)
        return all(self._parity(row & value) == 0 for row in self.h_rows)

    def decode_llr(self, llr: list[float]) -> tuple[list[int], list[int], bool, int]:
        if len(llr) != self.n:
            raise ValueError(f"LLR length must be {self.n}")
        q = [llr[i * self.z : (i + 1) * self.z] for i in range(self.nb)]
        r_msgs = {(i, j): [0.0] * self.z for i in range(self.mb) for j, _ in self.layers[i]}

        for it in range(1, self.max_iter + 1):
            for layer_idx, conns in enumerate(self.layers):
                extr: list[list[float]] = []
                for j, sh in conns:
                    v2c = [q[j][x] - r_msgs[(layer_idx, j)][x] for x in range(self.z)]
                    extr.append([v2c[(x + sh) % self.z] for x in range(self.z)])

                for lane in range(self.z):
                    signs = [-1.0 if extr[e][lane] < 0 else 1.0 for e in range(len(conns))]
                    sign_prod = 1.0
                    for s in signs:
                        sign_prod *= s
                    mags = [abs(extr[e][lane]) for e in range(len(conns))]
                    min1_idx = min(range(len(mags)), key=lambda e: mags[e])
                    min1 = mags[min1_idx]
                    min2 = min(mags[e] for e in range(len(mags)) if e != min1_idx)
                    for e, (j, sh) in enumerate(conns):
                        mag = min2 if e == min1_idx else min1
                        value = self.alpha * sign_prod * signs[e] * mag
                        dst = (lane + sh) % self.z
                        old = r_msgs[(layer_idx, j)][dst]
                        q[j][dst] += value - old
                        r_msgs[(layer_idx, j)][dst] = value

            hard = [1 if x < 0 else 0 for block in q for x in block]
            if self.check_bits(hard):
                return hard[: self.k], hard, True, it

        hard = [1 if x < 0 else 0 for block in q for x in block]
        return hard[: self.k], hard, self.check_bits(hard), self.max_iter

    def summary(self) -> str:
        return (
            f"QC-LDPC(N={self.n}, K={self.k}, M={self.m}, rate={self.k/self.n:.4f}, "
            f"Z={self.z}, PCM={self.mb}x{self.nb}, alpha={self.alpha}, max_iter={self.max_iter})"
        )


def pad_payload(payload: bytes, info_bytes: int) -> tuple[bytes, int]:
    pad_len = (-len(payload)) % info_bytes
    return payload + bytes(pad_len), pad_len


def encode_payload(payload: bytes, ldpc: QCLDPC, raw_blocks: bool) -> tuple[bytes, int]:
    if raw_blocks:
        if len(payload) % ldpc.info_bytes:
            raise ValueError(f"Encode input must be a multiple of {ldpc.info_bytes} bytes in raw block mode")
        data = payload
        pad_len = 0
    else:
        data, pad_len = pad_payload(payload, ldpc.info_bytes)
    out = bytearray()
    for i in range(0, len(data), ldpc.info_bytes):
        msg_bits = bytes_to_bits(data[i : i + ldpc.info_bytes])
        cw_bits = ldpc.encode_bits(msg_bits)
        out.extend(bits_to_bytes(cw_bits))
    return bytes(out), pad_len


def decode_payload(encoded: bytes, ldpc: QCLDPC, raw_blocks: bool, llr_mag: float) -> tuple[bytes, bool, int]:
    if len(encoded) % ldpc.codeword_bytes:
        raise ValueError(f"Decode input must be a multiple of {ldpc.codeword_bytes} bytes")
    decoded = bytearray()
    ok_all = True
    max_iter_used = 0
    for i in range(0, len(encoded), ldpc.codeword_bytes):
        hard_bits = bytes_to_bits(encoded[i : i + ldpc.codeword_bytes])
        llr = [llr_mag if b == 0 else -llr_mag for b in hard_bits]
        msg_bits, _, ok, it = ldpc.decode_llr(llr)
        decoded.extend(bits_to_bytes(msg_bits))
        ok_all = ok_all and ok
        max_iter_used = max(max_iter_used, it)
    return bytes(decoded), ok_all, max_iter_used


def inject_random_errors(encoded: bytes, codeword_bits: int, errors_per_codeword: int, seed: int) -> bytes:
    if errors_per_codeword <= 0:
        return encoded
    rng = random.Random(seed)
    bits = bytes_to_bits(encoded)
    for block_start in range(0, len(bits), codeword_bits):
        positions = rng.sample(range(codeword_bits), errors_per_codeword)
        for p in positions:
            bits[block_start + p] ^= 1
    return bits_to_bytes(bits)


def inject_manual_errors(encoded: bytes, codeword_bits: int, per_codeword_positions: list[int], global_positions: list[int]) -> bytes:
    if not per_codeword_positions and not global_positions:
        return encoded
    bits = bytes_to_bits(encoded)
    for pos in per_codeword_positions:
        if pos < 0 or pos >= codeword_bits:
            raise ValueError(f"Each --error-positions value must be in range 0..{codeword_bits-1}")
    for block_start in range(0, len(bits), codeword_bits):
        for pos in per_codeword_positions:
            bits[block_start + pos] ^= 1
    for pos in global_positions:
        if pos < 0 or pos >= len(bits):
            raise ValueError(f"Global error position {pos} is outside encoded bit range 0..{len(bits)-1}")
        bits[pos] ^= 1
    return bits_to_bytes(bits)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generic QC-LDPC hex encode/decode model")
    parser.add_argument("mode", choices=["encode", "decode", "roundtrip", "info"], help="Operation mode")
    parser.add_argument("input", nargs="?", type=Path, help="Input hex .txt file")
    parser.add_argument("--preset", choices=sorted(PRESETS), default="1120_840", help="Built-in LDPC preset")
    parser.add_argument("--config", type=Path, help="JSON config with z, pcm, alpha, max_iter")
    parser.add_argument("--z", type=int, help="Override Z")
    parser.add_argument("--alpha", type=float, help="Override NMSA alpha")
    parser.add_argument("--max-iter", type=int, help="Override decoder max iterations")
    parser.add_argument("--llr-mag", type=float, default=8.0, help="LLR magnitude used when decoding hard hex codewords")
    parser.add_argument("--enc-out", type=Path, default=Path("ldpc_encoded_hex.txt"), help="Encoded LDPC hex output")
    parser.add_argument("--dec-out", type=Path, default=Path("ldpc_decoded_hex.txt"), help="Decoded payload hex output")
    parser.add_argument("--raw-blocks", action="store_true", help="Do not add/remove length header or padding")
    parser.add_argument("--inject-errors", type=int, default=0, help="Random bit errors per LDPC codeword in roundtrip mode")
    parser.add_argument("--error-positions", default="", help="Manual bit positions to flip in every LDPC codeword")
    parser.add_argument("--global-error-positions", default="", help="Manual bit positions to flip over the whole encoded stream")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for error injection")
    args = parser.parse_args()

    cfg = load_config(args)
    ldpc = QCLDPC(cfg["pcm"], int(cfg["z"]), float(cfg.get("alpha", 0.75)), int(cfg.get("max_iter", 50)))

    if args.mode == "info":
        print(ldpc.summary())
        print(f"Info bytes/codeword : {ldpc.info_bytes}")
        print(f"Codeword bytes      : {ldpc.codeword_bytes}")
        return
    if args.input is None:
        raise SystemExit("input hex file is required for encode/decode/roundtrip")

    if args.mode == "encode":
        payload = read_hex_file(args.input)
        encoded, pad_len = encode_payload(payload, ldpc, args.raw_blocks)
        write_hex_file(args.enc_out, encoded)
        print("LDPC hex encode")
        print(ldpc.summary())
        print(f"Input hex file       : {args.input}")
        print(f"Output encoded file  : {args.enc_out}")
        print(f"Raw block mode       : {args.raw_blocks}")
        print(f"Payload bytes        : {len(payload)}")
        print(f"Zero padding bytes   : {pad_len}")
        print(f"Encoded bytes        : {len(encoded)}")
        print(f"Encoded hex          : {encoded.hex()}")

    elif args.mode == "decode":
        encoded = read_hex_file(args.input)
        decoded, ok, it = decode_payload(encoded, ldpc, args.raw_blocks, args.llr_mag)
        write_hex_file(args.dec_out, decoded)
        print("LDPC hex decode")
        print(ldpc.summary())
        print(f"Input encoded file   : {args.input}")
        print(f"Output decoded file  : {args.dec_out}")
        print(f"Raw block mode       : {args.raw_blocks}")
        print(f"Decode success       : {ok}")
        print(f"Max iter used        : {it}")
        print(f"Recovered hex        : {decoded.hex()}")

    else:
        payload = read_hex_file(args.input)
        encoded, pad_len = encode_payload(payload, ldpc, args.raw_blocks)
        encoded = inject_random_errors(encoded, ldpc.n, args.inject_errors, args.seed)
        manual_positions = parse_positions(args.error_positions)
        global_positions = parse_positions(args.global_error_positions)
        encoded = inject_manual_errors(encoded, ldpc.n, manual_positions, global_positions)
        write_hex_file(args.enc_out, encoded)
        decoded, ok, it = decode_payload(encoded, ldpc, args.raw_blocks, args.llr_mag)
        write_hex_file(args.dec_out, decoded)
        print("LDPC hex roundtrip")
        print(ldpc.summary())
        print(f"Input hex file       : {args.input}")
        print(f"Output encoded file  : {args.enc_out}")
        print(f"Output decoded file  : {args.dec_out}")
        print(f"Raw block mode       : {args.raw_blocks}")
        print(f"Zero padding bytes   : {pad_len}")
        print(f"Injected errors/cw   : {args.inject_errors}")
        print(f"Manual errors/cw     : {manual_positions}")
        print(f"Global manual errors : {global_positions}")
        print(f"Decode success       : {ok}")
        print(f"Max iter used        : {it}")
        print(f"Original hex         : {payload.hex()}")
        print(f"Recovered hex        : {decoded.hex()}")
        recovered_prefix_ok = decoded[:len(payload)] == payload
        padding_ok = all(x == 0 for x in decoded[len(payload):])
        print(f"Recovered prefix OK  : {recovered_prefix_ok}")
        print(f"Padding zero OK      : {padding_ok}")
        print(f"Roundtrip OK         : {recovered_prefix_ok and padding_ok and ok}")


if __name__ == "__main__":
    main()

