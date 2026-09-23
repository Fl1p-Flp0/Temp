#!/usr/bin/env python3
"""Pure-Python RX benchmark for the RTL AES-BCH-LDPC receive chain.

Measured path: 1120 signed 6-bit LLRs -> fixed-point QC-LDPC -> BCH(t=4)
-> padding check -> AES-128 decrypt x6. Corpus construction is completed
before timing; file I/O, DMA, Linux, and FPGA queue time are excluded.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import statistics
import struct
import sys
import time
from pathlib import Path

from benchmark_tx_python import (
    AES128,
    DEFAULT_KEY,
    DEFAULT_PAYLOAD,
    PAYLOAD_BYTES,
    RtlMatchedTxEncoder,
    parse_hex,
)


CODEWORD_BITS = 1120
BCH_BITS = 840
LLR_MAG_CLEAN = 12
HIL_PAYLOAD = bytes.fromhex(
    "00000000000000000000000000000000"
    "00112233445566778899aabbccddeeff"
    "0123456789abcdeffedcba9876543210"
    "112233445566778899aabbccddeeff00"
    "ffeeddccbbaa99887766554433221100"
    "13579bdf2468ace055aa55aaaa55aa55"
)
PCM = (
    (5, 14, 12, 1, 2, 37, 45, 26, 24, 0, 3, -1, 34, 7, 46, 10, -1, -1, -1, -1),
    (0, 35, 1, 26, 0, 10, 16, 16, 34, 4, 2, 23, 0, 51, -1, 49, 20, -1, -1, -1),
    (12, 28, 22, 46, 3, 16, 51, 2, 25, 29, 19, 18, 52, -1, 37, -1, 34, 39, -1, -1),
    (0, 51, 16, 31, 13, 39, 27, 33, 8, 27, 53, 13, -1, 52, 33, -1, -1, 38, 7, -1),
    (36, 6, 3, 51, 4, 19, 4, 45, 48, 9, -1, 11, 22, 23, 43, -1, -1, -1, 14, 1),
)


def sat(value: int) -> int:
    return max(-32, min(31, value))


def bits_from_bytes(data: bytes) -> list[int]:
    return [(byte >> shift) & 1 for byte in data for shift in range(7, -1, -1)]


def bytes_from_bits(bits: list[int]) -> bytes:
    if len(bits) % 8:
        raise ValueError("bit count must be byte aligned")
    return bytes(sum((bit & 1) << (7 - offset) for offset, bit in enumerate(bits[index : index + 8]))
                 for index in range(0, len(bits), 8))


class AES128Decrypt:
    """AES decrypt using the same key schedule and byte ordering as the RTL."""

    def __init__(self, key: bytes) -> None:
        encrypt = AES128(key)
        self.round_keys = encrypt.round_keys
        inverse = [0] * 256
        from benchmark_tx_python import SBOX
        for index, value in enumerate(SBOX):
            inverse[value] = index
        self.inv_sbox = tuple(inverse)

    @staticmethod
    def _inv_shift_rows(state: list[int]) -> list[int]:
        return [state[0], state[13], state[10], state[7], state[4], state[1], state[14], state[11],
                state[8], state[5], state[2], state[15], state[12], state[9], state[6], state[3]]

    @staticmethod
    def _mul(value: int, coefficient: int) -> int:
        result = 0
        while coefficient:
            if coefficient & 1:
                result ^= value
            value = ((value << 1) ^ (0x11B if value & 0x80 else 0)) & 0xFF
            coefficient >>= 1
        return result

    @classmethod
    def _inv_mix_columns(cls, state: list[int]) -> list[int]:
        mixed = [0] * 16
        for column in range(0, 16, 4):
            s0, s1, s2, s3 = state[column : column + 4]
            mixed[column] = cls._mul(s0, 14) ^ cls._mul(s1, 11) ^ cls._mul(s2, 13) ^ cls._mul(s3, 9)
            mixed[column + 1] = cls._mul(s0, 9) ^ cls._mul(s1, 14) ^ cls._mul(s2, 11) ^ cls._mul(s3, 13)
            mixed[column + 2] = cls._mul(s0, 13) ^ cls._mul(s1, 9) ^ cls._mul(s2, 14) ^ cls._mul(s3, 11)
            mixed[column + 3] = cls._mul(s0, 11) ^ cls._mul(s1, 13) ^ cls._mul(s2, 9) ^ cls._mul(s3, 14)
        return mixed

    def decrypt_block(self, ciphertext: bytes) -> bytes:
        if len(ciphertext) != 16:
            raise ValueError("AES ciphertext block must be 16 bytes")
        state = [left ^ right for left, right in zip(ciphertext, self.round_keys[10])]
        for round_index in range(9, 0, -1):
            state = [self.inv_sbox[value] for value in self._inv_shift_rows(state)]
            state = [left ^ right for left, right in zip(state, self.round_keys[round_index])]
            state = self._inv_mix_columns(state)
        state = [self.inv_sbox[value] for value in self._inv_shift_rows(state)]
        return bytes(left ^ right for left, right in zip(state, self.round_keys[0]))


class FixedPointQcLdpc:
    """Layered normalized Min-Sum matching the RTL's 6-bit arithmetic."""

    z = 56
    max_iter = 10

    def __init__(self) -> None:
        self.layers: list[list[tuple[int, int, int]]] = []
        edge_index = 0
        for row in PCM:
            layer: list[tuple[int, int, int]] = []
            for column, shift in enumerate(row):
                if shift >= 0:
                    layer.append((edge_index, column, shift))
                    edge_index += 1
            self.layers.append(layer)
        if edge_index != 79:
            raise RuntimeError("Unexpected QC-LDPC edge count")

    def syndrome_zero(self, q: list[int]) -> bool:
        for layer in self.layers:
            for z_index in range(self.z):
                parity = 0
                for _, column, shift in layer:
                    parity ^= int(q[column * self.z + ((z_index + shift) % self.z)] < 0)
                if parity:
                    return False
        return True

    def decode(self, llrs: list[int]) -> tuple[list[int], bool, int]:
        if len(llrs) != CODEWORD_BITS or any(value < -32 or value > 31 for value in llrs):
            raise ValueError("RX input must be exactly 1120 signed 6-bit LLR values")
        q = llrs[:]
        r_messages = [[0] * self.z for _ in range(79)]
        for iteration in range(1, self.max_iter + 1):
            for layer in self.layers:
                for z_index in range(self.z):
                    values: list[tuple[int, int, int]] = []
                    sign_all = 0
                    min1 = 63
                    min2 = 63
                    min1_pos = -1
                    for position, (edge, column, shift) in enumerate(layer):
                        lane = column * self.z + ((z_index + shift) % self.z)
                        v2c = sat(q[lane] - r_messages[edge][z_index])
                        magnitude = abs(v2c)
                        values.append((edge, lane, v2c))
                        sign_all ^= int(v2c < 0)
                        if magnitude < min1:
                            min2, min1, min1_pos = min1, magnitude, position
                        elif magnitude < min2:
                            min2 = magnitude
                    for position, (edge, lane, v2c) in enumerate(values):
                        magnitude = min2 if position == min1_pos else min1
                        scaled = min(31, (magnitude * 3) >> 2)
                        r_new = -scaled if (sign_all ^ int(v2c < 0)) else scaled
                        q[lane] = sat(v2c + r_new)
                        r_messages[edge][z_index] = r_new
            if self.syndrome_zero(q):
                return [int(value < 0) for value in q], True, iteration
        return [int(value < 0) for value in q], self.syndrome_zero(q), self.max_iter


class BchT4:
    """BCH(840,800,t=4) decoder over the RTL GF(2^10), x^10+x^7+1."""

    n_full = 1023
    t = 4

    def __init__(self) -> None:
        self.exp = [0] * (2 * self.n_full)
        self.log = [-1] * 1024
        value = 1
        for exponent in range(self.n_full):
            self.exp[exponent] = value
            self.log[value] = exponent
            value = self.mul(value, 2)
        if value != 1:
            raise RuntimeError("GF(2^10) primitive element check failed")
        for exponent in range(self.n_full, 2 * self.n_full):
            self.exp[exponent] = self.exp[exponent - self.n_full]

    @staticmethod
    def mul(left: int, right: int) -> int:
        product = 0
        for bit in range(10):
            if (right >> bit) & 1:
                product ^= left << bit
        for bit in range(18, 9, -1):
            if (product >> bit) & 1:
                product ^= 1 << bit
                product ^= 1 << (bit - 7)
                product ^= 1 << (bit - 10)
        return product & 0x3FF

    def divide(self, left: int, right: int) -> int:
        if left == 0:
            return 0
        if right == 0:
            raise ZeroDivisionError("GF division by zero")
        return self.exp[(self.log[left] - self.log[right]) % self.n_full]

    def syndromes(self, word: int) -> list[int]:
        values = [0] * (2 * self.t)
        for position in range(BCH_BITS):
            if (word >> position) & 1:
                for syndrome_index in range(1, 2 * self.t + 1):
                    values[syndrome_index - 1] ^= self.exp[(syndrome_index * position) % self.n_full]
        return values

    def berlekamp_massey(self, syndromes: list[int]) -> tuple[list[int], int]:
        sigma = [1] + [0] * (2 * self.t)
        b_poly = [1] + [0] * (2 * self.t)
        degree = 0
        shift = 1
        last_discrepancy = 1
        for index in range(2 * self.t):
            discrepancy = syndromes[index]
            for coefficient in range(1, degree + 1):
                discrepancy ^= self.mul(sigma[coefficient], syndromes[index - coefficient])
            if discrepancy == 0:
                shift += 1
                continue
            old_sigma = sigma[:]
            factor = self.divide(discrepancy, last_discrepancy)
            for coefficient in range(0, 2 * self.t + 1 - shift):
                if b_poly[coefficient]:
                    sigma[coefficient + shift] ^= self.mul(factor, b_poly[coefficient])
            if 2 * degree <= index:
                degree = index + 1 - degree
                b_poly = old_sigma
                last_discrepancy = discrepancy
                shift = 1
            else:
                shift += 1
        return sigma[: degree + 1], degree

    def decode(self, codeword: bytes) -> tuple[bytes, bytes, bool, int]:
        if len(codeword) != BCH_BITS // 8:
            raise ValueError("BCH input must be exactly 840 bits")
        word = int.from_bytes(codeword, "big")
        syndromes = self.syndromes(word)
        if not any(syndromes):
            return codeword[:100], codeword, True, 0
        sigma, degree = self.berlekamp_massey(syndromes)
        error_positions: list[int] = []
        for position in range(BCH_BITS):
            x_inverse = self.exp[(self.n_full - position) % self.n_full]
            value = 0
            power = 1
            for coefficient in sigma:
                if coefficient:
                    value ^= self.mul(coefficient, power)
                power = self.mul(power, x_inverse)
            if value == 0:
                error_positions.append(position)
        ok = degree <= self.t and len(error_positions) == degree
        corrected = word
        if ok:
            for position in error_positions:
                corrected ^= 1 << position
            ok = not any(self.syndromes(corrected))
        corrected_bytes = corrected.to_bytes(BCH_BITS // 8, "big")
        return corrected_bytes[:100], corrected_bytes, ok, len(error_positions) if ok else 0


class XorShift64Star:
    """The same deterministic RNG sequence used by the Linux HIL client."""

    mask = (1 << 64) - 1

    def __init__(self, seed: int) -> None:
        self.state = seed & self.mask
        if self.state == 0:
            self.state = 1

    def next(self) -> int:
        value = self.state
        value ^= value >> 12
        value ^= (value << 25) & self.mask
        value ^= value >> 27
        self.state = value & self.mask
        return (self.state * 2685821657736338717) & self.mask

    def uniform_open(self) -> float:
        return ((self.next() >> 11) + 0.5) / 9007199254740992.0

    def normal(self) -> float:
        return math.sqrt(-2.0 * math.log(self.uniform_open())) * math.cos(2.0 * math.pi * self.uniform_open())


def make_awgn_llrs(codeword: bytes, ebn0_db: float, llr_scale: float, rng: XorShift64Star) -> list[int]:
    rate = 768.0 / CODEWORD_BITS
    sigma = math.sqrt(1.0 / (2.0 * rate * (10.0 ** (ebn0_db / 10.0))))
    llrs: list[int] = []
    for bit in bits_from_bytes(codeword):
        sample = (-1.0 if bit else 1.0) + sigma * rng.normal()
        llrs.append(sat(int(round(llr_scale * (2.0 * sample / (sigma * sigma))))))
    return llrs


def clean_llrs(codeword: bytes) -> list[int]:
    return [-LLR_MAG_CLEAN if bit else LLR_MAG_CLEAN for bit in bits_from_bytes(codeword)]


def save_corpus(path: Path, corpus: list[list[int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    flattened = [value for llrs in corpus for value in llrs]
    path.write_bytes(struct.pack(f"{len(flattened)}b", *flattened))


def load_corpus(path: Path, frames: int) -> list[list[int]]:
    raw = path.read_bytes()
    expected = frames * CODEWORD_BITS
    if len(raw) != expected:
        raise ValueError(f"LLR corpus has {len(raw)} bytes; expected {expected} for {frames} frames")
    values = struct.unpack(f"{expected}b", raw)
    return [list(values[index : index + CODEWORD_BITS]) for index in range(0, expected, CODEWORD_BITS)]


def decode_frame(ldpc: FixedPointQcLdpc, bch: BchT4, aes: AES128Decrypt, llrs: list[int]) -> tuple[bytes, bool, int]:
    hard_bits, _syndrome_ok, iterations = ldpc.decode(llrs)
    bch_data, _corrected, bch_ok, _error_count = bch.decode(bytes_from_bits(hard_bits[:BCH_BITS]))
    # BCH is always attempted after LDPC reaches syndrome=0 or Imax.
    if not bch_ok or bch_data[96:] != bytes(4):
        return bytes(PAYLOAD_BYTES), False, iterations
    plaintext = b"".join(aes.decrypt_block(bch_data[index : index + 16]) for index in range(0, PAYLOAD_BYTES, 16))
    return plaintext, True, iterations


def validate_components(tx: RtlMatchedTxEncoder, ldpc: FixedPointQcLdpc, bch: BchT4, aes: AES128Decrypt) -> None:
    encrypted = AES128(DEFAULT_KEY).encrypt_block(DEFAULT_PAYLOAD[:16])
    if aes.decrypt_block(encrypted) != DEFAULT_PAYLOAD[:16]:
        raise RuntimeError("AES decrypt self-test failed")
    aes_frame = b"".join(tx.aes.encrypt_block(DEFAULT_PAYLOAD[index : index + 16]) for index in range(0, 96, 16))
    bch_word = tx.bch_encode(aes_frame)
    bch_word_int = int.from_bytes(bch_word, "big")
    for positions in ((839,), (839, 640), (839, 640, 320), (839, 640, 320, 0)):
        received = bch_word_int
        for position in positions:
            received ^= 1 << position
        data, corrected, ok, count = bch.decode(received.to_bytes(105, "big"))
        if not ok or count != len(positions) or corrected != bch_word or data != bch_word[:100]:
            raise RuntimeError(f"BCH t=4 self-test failed for {len(positions)} injected errors")
    received = bch_word_int
    for position in (839, 640, 320, 32, 0):
        received ^= 1 << position
    if bch.decode(received.to_bytes(105, "big"))[2]:
        raise RuntimeError("BCH t=4 self-test did not reject five injected errors")
    codeword = tx.encode_frame(DEFAULT_PAYLOAD)
    plaintext, accepted, _iterations = decode_frame(ldpc, bch, aes, clean_llrs(codeword))
    if not accepted or plaintext != DEFAULT_PAYLOAD:
        raise RuntimeError("Clean AES-BCH-LDPC RX self-test failed")


def write_csv(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def write_rows_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("At least one timed repetition is required")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark RX on a host CPU using a 6-bit fixed-point model")
    parser.add_argument("--frames", type=int, default=100, help="Measured RX frames")
    parser.add_argument("--warmup", type=int, default=5, help="Untimed RX frames before every timed repetition")
    parser.add_argument("--repetitions", type=int, default=7, help="Independent timed repetitions; median is the primary result")
    parser.add_argument("--scenario", choices=("clean", "awgn"), default="clean")
    parser.add_argument("--ebn0-db", type=float, default=4.0, help="AWGN Eb/N0 for --scenario awgn")
    parser.add_argument("--llr-scale", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20260922)
    parser.add_argument("--key-hex", default=DEFAULT_KEY.hex())
    parser.add_argument("--payload-hex", default=HIL_PAYLOAD.hex(), help="96-byte plaintext; default matches the board HIL client")
    parser.add_argument("--corpus", type=Path, help="Read exactly frames*1120 signed-LLR bytes from this file")
    parser.add_argument("--write-corpus", type=Path, help="Save the generated signed-LLR corpus for RV32/FPGA replay")
    parser.add_argument("--out-dir", type=Path, default=Path("build") / "cpu_benchmark" / "rx_python")
    parser.add_argument("--no-self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.frames <= 0 or args.warmup < 0 or args.repetitions <= 0 or args.llr_scale <= 0.0:
        raise SystemExit("frames and repetitions must be positive; warmup non-negative; llr-scale positive")
    key = parse_hex(args.key_hex, 16, "--key-hex")
    base_payload = parse_hex(args.payload_hex, PAYLOAD_BYTES, "--payload-hex")
    tx, ldpc, bch, aes = RtlMatchedTxEncoder(key), FixedPointQcLdpc(), BchT4(), AES128Decrypt(key)
    if not args.no_self_test:
        validate_components(tx, ldpc, bch, aes)

    # Keep one board-compatible frame across the corpus so every platform can
    # replay the identical LLR vectors and expected plaintext.
    expected_payloads = [base_payload] * args.frames
    if args.corpus:
        corpus = load_corpus(args.corpus, args.frames)
        corpus_source = str(args.corpus.resolve())
    else:
        rng = XorShift64Star(args.seed)
        corpus = []
        for payload in expected_payloads:
            codeword = tx.encode_frame(payload)
            corpus.append(clean_llrs(codeword) if args.scenario == "clean" else make_awgn_llrs(codeword, args.ebn0_db, args.llr_scale, rng))
        corpus_source = "generated"
    if args.write_corpus:
        save_corpus(args.write_corpus, corpus)

    repetition_rows: list[dict[str, object]] = []
    latency_values: list[float] = []
    reference_outcome: tuple[int, int, int, int, str] | None = None
    for repetition in range(1, args.repetitions + 1):
        for frame_index in range(args.warmup):
            decode_frame(ldpc, bch, aes, corpus[frame_index % args.frames])

        accepted = rejected = wrong_accept = total_iterations = 0
        output_hash = hashlib.sha256()
        total_start = time.perf_counter_ns()
        for expected, llrs in zip(expected_payloads, corpus):
            plaintext, ok, iterations = decode_frame(ldpc, bch, aes, llrs)
            total_iterations += iterations
            if not ok:
                rejected += 1
                output_hash.update(b"REJECT" + plaintext)
            elif plaintext != expected:
                wrong_accept += 1
                output_hash.update(b"WRONG!" + plaintext)
            else:
                accepted += 1
                output_hash.update(b"ACCEPT" + plaintext)
        elapsed_ns = time.perf_counter_ns() - total_start
        latency_us = elapsed_ns / args.frames / 1_000
        output_sha256 = output_hash.hexdigest()
        outcome = (accepted, rejected, wrong_accept, total_iterations, output_sha256)
        if reference_outcome is None:
            reference_outcome = outcome
        elif outcome != reference_outcome:
            raise RuntimeError("RX output or LDPC iteration count differs between repetitions")
        latency_values.append(latency_us)
        repetition_rows.append({
            "repetition": repetition,
            "frames": args.frames,
            "total_elapsed_s": f"{elapsed_ns / 1_000_000_000:.9f}",
            "rx_latency_avg_us_per_frame": f"{latency_us:.6f}",
            "ldpc_avg_iterations": f"{total_iterations / args.frames:.6f}",
            "accepted_frames": accepted,
            "rejected_frames": rejected,
            "wrong_accepted_frames": wrong_accept,
            "output_sha256": output_sha256,
        })

    if reference_outcome is None:
        raise RuntimeError("RX benchmark did not complete a timed repetition")
    accepted, rejected, wrong_accept, total_iterations, output_sha256 = reference_outcome

    row: dict[str, object] = {
        "benchmark": "rx_python_fixed6_ldpc1120_840_bch840_800_t4_aes128",
        "implementation": "pure_python_fixed_point_rtl_matched",
        "scenario": "external_corpus" if args.corpus else args.scenario,
        "ebn0_db": args.ebn0_db if not args.corpus else "external",
        "llr_scale": args.llr_scale if not args.corpus else "external",
        "seed": args.seed if not args.corpus else "external",
        "corpus": corpus_source,
        "payload_bits_per_frame": 768,
        "llr_values_per_frame": CODEWORD_BITS,
        "frames_per_repetition": args.frames,
        "warmup_frames_per_repetition": args.warmup,
        "repetitions": args.repetitions,
        "primary_latency_statistic": "median_of_repetition_averages",
        "rx_latency_median_us_per_frame": f"{statistics.median(latency_values):.6f}",
        "rx_latency_mean_us_per_frame": f"{statistics.fmean(latency_values):.6f}",
        "rx_latency_stddev_us_per_frame": f"{statistics.pstdev(latency_values):.6f}",
        "rx_latency_min_us_per_frame": f"{min(latency_values):.6f}",
        "rx_latency_max_us_per_frame": f"{max(latency_values):.6f}",
        "ldpc_avg_iterations": f"{total_iterations / args.frames:.6f}",
        "accepted_frames_per_repetition": accepted,
        "rejected_frames_per_repetition": rejected,
        "wrong_accepted_frames_per_repetition": wrong_accept,
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "system": platform.platform(),
        "output_sha256": output_sha256,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "rx_python_summary.csv"
    runs_csv_path = args.out_dir / "rx_python_runs.csv"
    json_path = args.out_dir / "rx_python_summary.json"
    write_csv(csv_path, row)
    write_rows_csv(runs_csv_path, repetition_rows)
    json_path.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")
    print("PASS: CPU RX benchmark completed")
    print(f"scenario={row['scenario']} repetitions={args.repetitions} frames_per_repetition={args.frames} accepted={accepted} rejected={rejected} wrong_accept={wrong_accept}")
    print(f"RX primary latency (median) = {float(row['rx_latency_median_us_per_frame']):.3f} us/frame; LDPC average iterations = {float(row['ldpc_avg_iterations']):.3f}")
    print(f"RX mean +/- population stddev = {float(row['rx_latency_mean_us_per_frame']):.3f} +/- {float(row['rx_latency_stddev_us_per_frame']):.3f} us/frame")
    print(f"CSV={csv_path.resolve()}")
    print(f"RUNS_CSV={runs_csv_path.resolve()}")
    print(f"JSON={json_path.resolve()}")
    if args.write_corpus:
        print(f"CORPUS={args.write_corpus.resolve()}")


if __name__ == "__main__":
    main()
