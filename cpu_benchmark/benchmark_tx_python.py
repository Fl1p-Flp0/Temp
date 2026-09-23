#!/usr/bin/env python3
"""Pure-Python TX benchmark for AES-128 + BCH(840,800,t=4) + QC-LDPC(1120,840).

Measured path: 96-byte payload -> AES-128 x6 -> BCH -> LDPC -> 1120-bit word.
File I/O, AXI DMA, Linux driver, and FPGA queue time are excluded. BCH and
LDPC masks are embedded from the verified RTL, so this script is standalone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import statistics
import sys
import time
from pathlib import Path

from rtl_masks import BCH_PARITY_MASKS, LDPC_PARITY_MASKS


PAYLOAD_BYTES = 96
AES_BLOCK_BYTES = 16
BCH_CODEWORD_BITS = 840
LDPC_CODEWORD_BITS = 1120
DEFAULT_KEY = bytes.fromhex("000102030405060708090a0b0c0d0e0f")
DEFAULT_PAYLOAD = bytes.fromhex(
    "00112233445566778899aabbccddeeff"
    "102132435465768798a9bacbdcedfe0f"
    "0f1e2d3c4b5a69788796a5b4c3d2e1f0"
    "ffeeddccbbaa99887766554433221100"
    "0123456789abcdeffedcba9876543210"
    "55aa55aa0f0ff0f0a5a55a5ac3c33c3c"
)


def xtime(value: int) -> int:
    value <<= 1
    return (value ^ 0x11B) & 0xFF if value & 0x100 else value


def gf_mul(left: int, right: int) -> int:
    product = 0
    while right:
        if right & 1:
            product ^= left
        left = xtime(left)
        right >>= 1
    return product


def gf_pow(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = gf_mul(result, value)
        value = gf_mul(value, value)
        exponent >>= 1
    return result


def make_sbox() -> tuple[int, ...]:
    table: list[int] = []
    for value in range(256):
        inverse = 0 if value == 0 else gf_pow(value, 254)
        transformed = inverse
        for shift in range(1, 5):
            transformed ^= ((inverse << shift) | (inverse >> (8 - shift))) & 0xFF
        table.append((transformed ^ 0x63) & 0xFF)
    return tuple(table)


SBOX = make_sbox()
RCON = (0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1B, 0x36)


class AES128:
    """FIPS-197 byte order; matches the AES RTL core."""

    def __init__(self, key: bytes) -> None:
        if len(key) != AES_BLOCK_BYTES:
            raise ValueError("AES-128 key must be exactly 16 bytes")
        self.round_keys = self._expand_key(key)

    @staticmethod
    def _expand_key(key: bytes) -> tuple[bytes, ...]:
        expanded = list(key)
        word_index = 4
        while len(expanded) < 176:
            temp = expanded[-4:]
            if word_index % 4 == 0:
                temp = [SBOX[temp[1]], SBOX[temp[2]], SBOX[temp[3]], SBOX[temp[0]]]
                temp[0] ^= RCON[word_index // 4]
            for offset in range(4):
                expanded.append(expanded[-16] ^ temp[offset])
            word_index += 1
        return tuple(bytes(expanded[index : index + 16]) for index in range(0, 176, 16))

    @staticmethod
    def _shift_rows(state: list[int]) -> list[int]:
        return [state[0], state[5], state[10], state[15], state[4], state[9], state[14], state[3],
                state[8], state[13], state[2], state[7], state[12], state[1], state[6], state[11]]

    @staticmethod
    def _mix_columns(state: list[int]) -> list[int]:
        mixed = [0] * 16
        for column in range(0, 16, 4):
            s0, s1, s2, s3 = state[column : column + 4]
            m20, m21, m22, m23 = xtime(s0), xtime(s1), xtime(s2), xtime(s3)
            mixed[column] = m20 ^ (m21 ^ s1) ^ s2 ^ s3
            mixed[column + 1] = s0 ^ m21 ^ (m22 ^ s2) ^ s3
            mixed[column + 2] = s0 ^ s1 ^ m22 ^ (m23 ^ s3)
            mixed[column + 3] = (m20 ^ s0) ^ s1 ^ s2 ^ m23
        return mixed

    def encrypt_block(self, plaintext: bytes) -> bytes:
        if len(plaintext) != AES_BLOCK_BYTES:
            raise ValueError("AES plaintext block must be exactly 16 bytes")
        state = [left ^ right for left, right in zip(plaintext, self.round_keys[0])]
        for round_index in range(1, 11):
            state = self._shift_rows([SBOX[value] for value in state])
            if round_index != 10:
                state = self._mix_columns(state)
            state = [left ^ right for left, right in zip(state, self.round_keys[round_index])]
        return bytes(state)


class RtlMatchedTxEncoder:
    """Combinational model of the BCH and LDPC TX RTL datapaths."""

    def __init__(self, key: bytes) -> None:
        self.aes = AES128(key)
        self.bch_masks = BCH_PARITY_MASKS
        self.ldpc_masks = LDPC_PARITY_MASKS

    def bch_encode(self, aes_frame: bytes) -> bytes:
        if len(aes_frame) != PAYLOAD_BYTES:
            raise ValueError("AES frame must contain six blocks")
        data_value = int.from_bytes(aes_frame + bytes(4), "big")  # {aes_frame, 32'b0}
        parity = 0
        remaining = data_value
        while remaining:
            lowest = remaining & -remaining
            parity ^= self.bch_masks[lowest.bit_length() - 1]
            remaining ^= lowest
        return ((data_value << 40) | parity).to_bytes(BCH_CODEWORD_BITS // 8, "big")

    def ldpc_encode(self, bch_codeword: bytes) -> bytes:
        if len(bch_codeword) != BCH_CODEWORD_BITS // 8:
            raise ValueError("BCH codeword must contain 840 bits")
        data_value = int.from_bytes(bch_codeword, "big")
        parity = 0
        for mask in self.ldpc_masks:
            parity = (parity << 1) | ((data_value & mask).bit_count() & 1)
        return ((data_value << 280) | parity).to_bytes(LDPC_CODEWORD_BITS // 8, "big")

    def encode_frame(self, plaintext: bytes) -> bytes:
        if len(plaintext) != PAYLOAD_BYTES:
            raise ValueError("TX payload must be exactly 96 bytes / 768 bits")
        aes_frame = b"".join(self.aes.encrypt_block(plaintext[index : index + 16]) for index in range(0, 96, 16))
        return self.ldpc_encode(self.bch_encode(aes_frame))


def payload_for_frame(base: bytes, frame_index: int) -> bytes:
    """Hold a 96-byte contract while avoiding a constant input stream."""
    suffix = frame_index.to_bytes(8, "big", signed=False)
    return base[:-8] + bytes(left ^ right for left, right in zip(base[-8:], suffix))


def parse_hex(text: str, expected_bytes: int, name: str) -> bytes:
    clean = "".join(text.split()).removeprefix("0x").removeprefix("0X")
    try:
        value = bytes.fromhex(clean)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"{name} must be hexadecimal") from exc
    if len(value) != expected_bytes:
        raise argparse.ArgumentTypeError(f"{name} must contain exactly {expected_bytes} bytes")
    return value


def validate_encoder(encoder: RtlMatchedTxEncoder) -> None:
    aes_known = bytes.fromhex("69c4e0d86a7b0430d8cdb78070b4c55a")
    aes_input = bytes.fromhex("00112233445566778899aabbccddeeff")
    if encoder.aes.encrypt_block(aes_input) != aes_known:
        raise RuntimeError("AES-128 self-test failed")
    codeword = encoder.encode_frame(DEFAULT_PAYLOAD)
    if len(codeword) != LDPC_CODEWORD_BITS // 8:
        raise RuntimeError("TX output width self-test failed")
    aes_frame = b"".join(encoder.aes.encrypt_block(DEFAULT_PAYLOAD[i : i + 16]) for i in range(0, 96, 16))
    if codeword[: BCH_CODEWORD_BITS // 8] != encoder.bch_encode(aes_frame):
        raise RuntimeError("LDPC systematic-prefix self-test failed")


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
    parser = argparse.ArgumentParser(description="Benchmark TX on a host CPU using a pure-Python RTL-matched model")
    parser.add_argument("--frames", type=int, default=10_000, help="Measured frames; each carries 768 useful bits")
    parser.add_argument("--warmup", type=int, default=250, help="Untimed warm-up frames before every timed repetition")
    parser.add_argument("--repetitions", type=int, default=7, help="Independent timed repetitions; median is the primary result")
    parser.add_argument("--key-hex", default=DEFAULT_KEY.hex(), help="AES-128 key, 16 bytes in hexadecimal")
    parser.add_argument("--payload-hex", default=DEFAULT_PAYLOAD.hex(), help="One 96-byte payload in hexadecimal")
    parser.add_argument("--out-dir", type=Path, default=Path("build") / "cpu_benchmark" / "tx_python")
    parser.add_argument("--emit-codeword", type=Path, help="Write one deterministic 1120-bit TX codeword as hexadecimal")
    parser.add_argument("--no-self-test", action="store_true", help="Skip setup-time AES and framing checks")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.frames <= 0 or args.warmup < 0 or args.repetitions <= 0:
        raise SystemExit("frames and repetitions must be positive; warmup must be non-negative")

    encoder = RtlMatchedTxEncoder(parse_hex(args.key_hex, 16, "--key-hex"))
    base_payload = parse_hex(args.payload_hex, PAYLOAD_BYTES, "--payload-hex")
    if not args.no_self_test:
        validate_encoder(encoder)

    # Do not charge deterministic test-vector generation to the TX datapath.
    payloads = [payload_for_frame(base_payload, frame_index) for frame_index in range(args.frames)]
    repetition_rows: list[dict[str, object]] = []
    latency_values: list[float] = []
    output_hashes: list[str] = []
    for repetition in range(1, args.repetitions + 1):
        for frame_index in range(args.warmup):
            encoder.encode_frame(payload_for_frame(base_payload, frame_index))

        output_hash = hashlib.sha256()
        total_start = time.perf_counter_ns()
        for payload in payloads:
            codeword = encoder.encode_frame(payload)
            output_hash.update(codeword)
        total_elapsed_ns = time.perf_counter_ns() - total_start
        latency_us = total_elapsed_ns / args.frames / 1_000
        output_sha256 = output_hash.hexdigest()
        latency_values.append(latency_us)
        output_hashes.append(output_sha256)
        repetition_rows.append({
            "repetition": repetition,
            "frames": args.frames,
            "total_elapsed_s": f"{total_elapsed_ns / 1_000_000_000:.9f}",
            "tx_latency_avg_us_per_frame": f"{latency_us:.6f}",
            "output_sha256": output_sha256,
        })

    if len(set(output_hashes)) != 1:
        raise RuntimeError("TX output differs between repetitions")

    if args.emit_codeword:
        args.emit_codeword.parent.mkdir(parents=True, exist_ok=True)
        args.emit_codeword.write_text(encoder.encode_frame(base_payload).hex() + "\n", encoding="ascii")

    median_latency_us = statistics.median(latency_values)
    payload_throughput_mbit_s = 768.0 / median_latency_us
    row: dict[str, object] = {
        "benchmark": "tx_python_aes128_bch840_800_t4_ldpc1120_840",
        "implementation": "pure_python_rtl_matched_masks",
        "payload_bits_per_frame": 768,
        "codeword_bits_per_frame": 1120,
        "frames_per_repetition": args.frames,
        "warmup_frames_per_repetition": args.warmup,
        "repetitions": args.repetitions,
        "primary_latency_statistic": "median_of_repetition_averages",
        "tx_latency_median_us_per_frame": f"{median_latency_us:.6f}",
        "tx_latency_mean_us_per_frame": f"{statistics.fmean(latency_values):.6f}",
        "tx_latency_stddev_us_per_frame": f"{statistics.pstdev(latency_values):.6f}",
        "tx_latency_min_us_per_frame": f"{min(latency_values):.6f}",
        "tx_latency_max_us_per_frame": f"{max(latency_values):.6f}",
        "payload_throughput_mbit_s": f"{payload_throughput_mbit_s:.6f}",
        "python": sys.version.split()[0],
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown",
        "system": platform.platform(),
        "output_sha256": output_hashes[0],
        "bch_parity_masks": "embedded_from_verified_rtl",
        "ldpc_parity_masks": "embedded_from_verified_rtl",
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "tx_python_summary.csv"
    runs_csv_path = args.out_dir / "tx_python_runs.csv"
    json_path = args.out_dir / "tx_python_summary.json"
    write_csv(csv_path, row)
    write_rows_csv(runs_csv_path, repetition_rows)
    json_path.write_text(json.dumps(row, indent=2) + "\n", encoding="utf-8")

    print("PASS: CPU TX benchmark completed")
    print(f"repetitions={args.repetitions} frames_per_repetition={args.frames} payload_bits_per_frame=768 codeword_bits_per_frame=1120")
    print(f"TX primary latency (median) = {float(row['tx_latency_median_us_per_frame']):.3f} us/frame")
    print(f"TX payload throughput = {float(row['payload_throughput_mbit_s']):.3f} Mbit/s")
    print(f"TX mean +/- population stddev = {float(row['tx_latency_mean_us_per_frame']):.3f} +/- {float(row['tx_latency_stddev_us_per_frame']):.3f} us/frame")
    print(f"CSV={csv_path.resolve()}")
    print(f"RUNS_CSV={runs_csv_path.resolve()}")
    print(f"JSON={json_path.resolve()}")
    if args.emit_codeword:
        print(f"CODEWORD={args.emit_codeword.resolve()}")


if __name__ == "__main__":
    main()
