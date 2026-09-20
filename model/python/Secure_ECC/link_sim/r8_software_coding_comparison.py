#!/usr/bin/env python3
"""Paired software comparison for the R8 coding study.

Study modes:

1. uncoded
   640 payload bits -> BPSK/AWGN -> hard decision.

2. ldpc_1120_640
   640 payload bits -> experimental QC-LDPC(1120,640).

3. bch_ldpc_fallback
   Algorithm-level BCH(840,640) + LDPC(1120,840). BCH receives the LDPC
   hard output even when the LDPC parity check fails. This mode studies the
   theoretical outer-code cleanup effect; it is not the present RTL policy.

4. bch_ldpc_r8
   640 payload bits -> BCH(840,640) -> the implemented
   QC-LDPC(1120,840) -> R8-equivalent fail-closed receive policy.

For every Eb/N0, seed and frame, all schemes use the same payload.  The two
1120-bit coded schemes also use the same 1120 standard-normal noise samples.
This common-random-number design reduces comparison variance without claiming
that different codewords produce identical received LLRs.

Both BCH-LDPC paths quantize channel LLRs to signed six-bit integers and use a
layered integer Normalized Min-Sum decoder with alpha=3/4 and ten iterations.
Only bch_ldpc_r8 applies the RTL gate that blocks BCH unless the LDPC parity
check succeeds. Rejected frames contribute 640 conservative effective
output-bit errors, matching the fail-closed HIL reporting convention.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import struct
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

from secure_ecc_link_sim import BCHShortened, count_errors, import_ldpc_model, parse_snr_list, rand_bits


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_LDPC640_CONFIG = THIS_DIR / "configs" / "ldpc_1120_640_config.json"
THEORY_MODES = ("uncoded", "ldpc_1120_640", "bch_ldpc_fallback")
R8_MODES = ("uncoded", "ldpc_1120_640", "bch_ldpc_r8")
AVAILABLE_MODES = (
    "uncoded",
    "ldpc_1120_640",
    "ldpc_640pad",
    "bch_ldpc_fallback",
    "bch_ldpc_r8",
)

MODE_LABELS = {
    "uncoded": "Uncoded BPSK",
    "ldpc_1120_640": "Experimental QC-LDPC(1120,640)",
    "ldpc_640pad": "LDPC(1120,840), 640 payload + 200 zero padding",
    "bch_ldpc_fallback": "BCH(840,640) + LDPC(1120,840), BCH fallback model",
    "bch_ldpc_r8": "R8 BCH(840,640) + LDPC(1120,840)",
}


def sat6(value: int) -> int:
    return max(-32, min(31, int(value)))


def quantize_llr6(value: float) -> int:
    # Python round and C lrint both use nearest-even with the default rounding
    # mode.  Clip after rounding to match secc_queue_dma_awgn_e2e.c.
    return sat6(round(value))


def bits_to_bytes(bits: Iterable[int]) -> bytes:
    values = list(bits)
    if len(values) % 8:
        raise ValueError("Bit vector length must be divisible by eight")
    output = bytearray()
    for offset in range(0, len(values), 8):
        byte = 0
        for bit in values[offset : offset + 8]:
            byte = (byte << 1) | (int(bit) & 1)
        output.append(byte)
    return bytes(output)


def parse_int_list(text: str) -> list[int]:
    values = [int(item.strip(), 0) for item in text.split(",") if item.strip()]
    if not values:
        raise ValueError("At least one seed is required")
    return values


def make_rng(seed: int, ebn0_db: float) -> random.Random:
    # String seeding is deterministic across supported CPython versions and
    # avoids dependence on process-randomized hash().
    return random.Random(f"r8-comparison:{seed}:{ebn0_db:.9f}")


def channel_llr6(
    bits: list[int],
    ebn0_db: float,
    rate: float,
    standard_noise: list[float],
    llr_scale: float,
) -> tuple[list[int], int]:
    if len(bits) > len(standard_noise):
        raise ValueError("Noise vector is shorter than the transmitted frame")
    ebn0_linear = 10.0 ** (ebn0_db / 10.0)
    sigma = math.sqrt(1.0 / (2.0 * rate * ebn0_linear))
    inv_sigma2_times_2 = 2.0 / (sigma * sigma)
    llr: list[int] = []
    channel_errors = 0
    for index, bit in enumerate(bits):
        symbol = 1.0 if bit == 0 else -1.0
        sample = symbol + sigma * standard_noise[index]
        channel_hard = 0 if sample >= 0.0 else 1
        channel_errors += int(channel_hard != bit)
        llr.append(quantize_llr6(llr_scale * inv_sigma2_times_2 * sample))
    return llr, channel_errors


class R8IntegerLayeredDecoder:
    """Signed-six-bit layered NMS model matching the R8 RTL arithmetic."""

    def __init__(self, code: Any, max_iter: int = 10) -> None:
        if abs(float(code.alpha) - 0.75) > 1e-12:
            raise ValueError("The R8 integer decoder requires alpha=3/4")
        self.code = code
        self.max_iter = int(max_iter)
        self.n = int(code.n)
        self.k = int(code.k)
        self.z = int(code.z)
        self.layers: list[list[tuple[int, int, int]]] = []
        edge_index = 0
        for connections in code.layers:
            layer: list[tuple[int, int, int]] = []
            for column, shift in connections:
                layer.append((int(column), int(shift), edge_index))
                edge_index += 1
            self.layers.append(layer)
        self.active_edges = edge_index

    def decode(self, channel_llr: list[int]) -> tuple[list[int], list[int], bool, int, list[int]]:
        if len(channel_llr) != self.n:
            raise ValueError(f"Expected {self.n} LLR values, got {len(channel_llr)}")

        lq = [sat6(value) for value in channel_llr]
        r_messages = [[0] * self.z for _ in range(self.active_edges)]

        for iteration in range(1, self.max_iter + 1):
            for layer in self.layers:
                for z_index in range(self.z):
                    captured: list[tuple[int, int, int, int]] = []
                    sign_all = 0
                    min1 = 63
                    min2 = 63
                    min1_column = -1

                    for column, shift, edge_index in layer:
                        shifted_z = (z_index + shift) % self.z
                        lane_address = column * self.z + shifted_z
                        old_r = r_messages[edge_index][z_index]
                        v2c = sat6(lq[lane_address] - old_r)
                        captured.append((column, edge_index, lane_address, v2c))

                        sign_all ^= int(v2c < 0)
                        magnitude = abs(v2c)
                        if magnitude < min1:
                            min2 = min1
                            min1 = magnitude
                            min1_column = column
                        elif magnitude < min2:
                            min2 = magnitude

                    for column, edge_index, lane_address, v2c in captured:
                        magnitude = min2 if column == min1_column else min1
                        scaled = min(31, (3 * magnitude) >> 2)
                        output_sign = sign_all ^ int(v2c < 0)
                        new_r = -scaled if output_sign else scaled
                        r_messages[edge_index][z_index] = new_r
                        lq[lane_address] = sat6(v2c + new_r)

            hard = [1 if value < 0 else 0 for value in lq]
            if self.code.check_bits(hard):
                return hard[: self.k], hard, True, iteration, lq

        hard = [1 if value < 0 else 0 for value in lq]
        return hard[: self.k], hard, self.code.check_bits(hard), self.max_iter, lq


def make_ldpc_1120_840(max_iter: int):
    module = import_ldpc_model()
    config = dict(module.PRESETS["1120_840"])
    return module.QCLDPC(config["pcm"], int(config["z"]), 0.75, int(max_iter))


def make_ldpc_1120_640(path: Path, max_iter: int):
    module = import_ldpc_model()
    config = json.loads(path.read_text(encoding="utf-8-sig"))
    code = module.QCLDPC(config["pcm"], int(config["z"]), 0.75, int(max_iter))
    if code.n != 1120 or code.k != 640:
        raise ValueError(f"Expected LDPC(1120,640), got ({code.n},{code.k})")
    return code, str(config.get("name", path.stem)), str(config.get("description", ""))


@dataclass
class FrameContext:
    ebn0_db: float
    seed: int
    frame_index: int
    payload: list[int]
    standard_noise: list[float]


def base_result(context: FrameContext, mode: str) -> dict[str, Any]:
    return {
        "mode": mode,
        "ebn0_db": context.ebn0_db,
        "seed": context.seed,
        "frame_index": context.frame_index,
        "payload_bits": 640,
        "transmitted_bits": 0,
        "channel_bit_errors": 0,
        "payload_bit_errors_valid": 0,
        "effective_output_bit_errors": 0,
        "frame_error": 0,
        "output_valid": 1,
        "decode_reject": 0,
        "accepted_wrong": 0,
        "ldpc_ok": "",
        "bch_ran": 0,
        "bch_ok": "",
        "ldpc_iterations": 0,
        "ldpc_residual_errors_840": "",
        "residual_class": "",
        "before_bch_payload_errors": "",
        "after_bch_payload_errors": "",
        "after_bch_output_valid": "",
        "bch_corrected_errors": 0,
        "bch_cleaned_frame": 0,
    }


def simulate_uncoded(context: FrameContext) -> dict[str, Any]:
    result = base_result(context, "uncoded")
    rate = 1.0
    ebn0_linear = 10.0 ** (context.ebn0_db / 10.0)
    sigma = math.sqrt(1.0 / (2.0 * rate * ebn0_linear))
    hard: list[int] = []
    for index, bit in enumerate(context.payload):
        symbol = 1.0 if bit == 0 else -1.0
        sample = symbol + sigma * context.standard_noise[index]
        hard.append(0 if sample >= 0.0 else 1)
    errors = count_errors(context.payload, hard)
    result.update(
        transmitted_bits=640,
        channel_bit_errors=errors,
        payload_bit_errors_valid=errors,
        effective_output_bit_errors=errors,
        frame_error=int(errors > 0),
        accepted_wrong=int(errors > 0),
    )
    return result


def simulate_ldpc_only(
    context: FrameContext,
    code: Any,
    decoder: R8IntegerLayeredDecoder,
    llr_scale: float,
    mode: str,
    padded: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    result = base_result(context, mode)
    message = context.payload + ([0] * 200 if padded else [])
    codeword = code.encode_bits(message)
    llr, channel_errors = channel_llr6(
        codeword, context.ebn0_db, 640 / 1120, context.standard_noise, llr_scale
    )
    decoded, hard, ok, iterations, _lq = decoder.decode(llr)
    payload_out = decoded[:640]
    payload_errors = count_errors(context.payload, payload_out)
    valid = bool(ok)
    result.update(
        transmitted_bits=1120,
        channel_bit_errors=channel_errors,
        payload_bit_errors_valid=payload_errors if valid else 0,
        effective_output_bit_errors=payload_errors if valid else 640,
        frame_error=int((not valid) or payload_errors > 0),
        output_valid=int(valid),
        decode_reject=int(not valid),
        accepted_wrong=int(valid and payload_errors > 0),
        ldpc_ok=int(valid),
        ldpc_iterations=iterations,
    )
    artifacts = {
        "codeword": codeword,
        "llr": llr,
        "hard": hard,
        "payload_out": payload_out if valid else None,
    }
    return result, artifacts


def residual_class(errors: int) -> str:
    if errors == 0:
        return "zero"
    if errors <= 20:
        return "1_to_20"
    return "gt_20"


def simulate_bch_ldpc(
    context: FrameContext,
    ldpc: Any,
    decoder: R8IntegerLayeredDecoder,
    bch: BCHShortened,
    llr_scale: float,
    mode: str,
    require_ldpc_ok: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    result = base_result(context, mode)
    bch_codeword = bch.encode_bits(context.payload)
    ldpc_codeword = ldpc.encode_bits(bch_codeword)
    llr, channel_errors = channel_llr6(
        ldpc_codeword, context.ebn0_db, 640 / 1120, context.standard_noise, llr_scale
    )
    ldpc_decoded, hard, ldpc_ok, iterations, _lq = decoder.decode(llr)
    ldpc_message = ldpc_decoded[:840]
    residual_errors = count_errors(bch_codeword, ldpc_message)
    before_payload_errors = count_errors(context.payload, ldpc_message[:640])

    result.update(
        transmitted_bits=1120,
        channel_bit_errors=channel_errors,
        ldpc_ok=int(ldpc_ok),
        ldpc_iterations=iterations,
        ldpc_residual_errors_840=residual_errors,
        residual_class=residual_class(residual_errors),
        before_bch_payload_errors=before_payload_errors,
    )

    bch_output: list[int] | None = None
    bch_ok: bool | None = None
    if ldpc_ok or not require_ldpc_ok:
        bch_output, bch_ok, corrected = bch.decode_bits(ldpc_message)
        after_errors = count_errors(context.payload, bch_output)
        valid = bool(bch_ok)
        result.update(
            bch_ran=1,
            bch_ok=int(valid),
            bch_corrected_errors=corrected,
            after_bch_payload_errors=after_errors,
            after_bch_output_valid=int(valid),
            bch_cleaned_frame=int(before_payload_errors > 0 and after_errors == 0 and valid),
            payload_bit_errors_valid=after_errors if valid else 0,
            effective_output_bit_errors=after_errors if valid else 640,
            frame_error=int((not valid) or after_errors > 0),
            output_valid=int(valid),
            decode_reject=int(not valid),
            accepted_wrong=int(valid and after_errors > 0),
        )
    else:
        # Match secure_ecc_rx_ctrl: in the deployed R8 policy BCH is not
        # started unless the LDPC parity check succeeds for the current frame.
        result.update(
            bch_ran=0,
            bch_ok="",
            after_bch_payload_errors="",
            after_bch_output_valid="",
            effective_output_bit_errors=640,
            frame_error=1,
            output_valid=0,
            decode_reject=1,
        )

    artifacts = {
        "bch_codeword": bch_codeword,
        "ldpc_codeword": ldpc_codeword,
        "llr": llr,
        "hard": hard,
        "payload_out": bch_output if bool(bch_ok) else None,
    }
    return result, artifacts


def simulate_bch_ldpc_r8(
    context: FrameContext,
    ldpc: Any,
    decoder: R8IntegerLayeredDecoder,
    bch: BCHShortened,
    llr_scale: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    return simulate_bch_ldpc(
        context,
        ldpc,
        decoder,
        bch,
        llr_scale,
        mode="bch_ldpc_r8",
        require_ldpc_ok=True,
    )


def simulate_bch_ldpc_fallback(
    context: FrameContext,
    ldpc: Any,
    decoder: R8IntegerLayeredDecoder,
    bch: BCHShortened,
    llr_scale: float,
) -> tuple[dict[str, Any], dict[str, Any]]:
    return simulate_bch_ldpc(
        context,
        ldpc,
        decoder,
        bch,
        llr_scale,
        mode="bch_ldpc_fallback",
        require_ldpc_ok=False,
    )


def wilson_interval(events: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0:
        return 0.0, 0.0
    p = events / trials
    denominator = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt((p * (1.0 - p) + z * z / (4.0 * trials)) / trials) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def make_aggregator() -> dict[str, float]:
    return defaultdict(float)


def update_aggregator(aggregate: dict[str, float], row: dict[str, Any]) -> None:
    aggregate["frames"] += 1
    aggregate["payload_bits"] += int(row["payload_bits"])
    aggregate["transmitted_bits"] += int(row["transmitted_bits"])
    aggregate["channel_bit_errors"] += int(row["channel_bit_errors"])
    aggregate["payload_bit_errors_valid"] += int(row["payload_bit_errors_valid"])
    aggregate["effective_output_bit_errors"] += int(row["effective_output_bit_errors"])
    if row["after_bch_payload_errors"] != "":
        aggregate["after_bch_payload_errors"] += int(row["after_bch_payload_errors"])
        aggregate["after_bch_payload_bits"] += 640
    aggregate["frame_errors"] += int(row["frame_error"])
    aggregate["valid_frames"] += int(row["output_valid"])
    aggregate["decode_rejects"] += int(row["decode_reject"])
    aggregate["accepted_wrong_frames"] += int(row["accepted_wrong"])
    aggregate["ldpc_iterations"] += int(row["ldpc_iterations"])
    aggregate["bch_ran_frames"] += int(row["bch_ran"])
    aggregate["bch_corrected_errors"] += int(row["bch_corrected_errors"])
    aggregate["bch_cleaned_frames"] += int(row["bch_cleaned_frame"])
    if row["residual_class"]:
        aggregate[f"residual_{row['residual_class']}_frames"] += 1
    if row["before_bch_payload_errors"] != "":
        aggregate["before_bch_payload_errors"] += int(row["before_bch_payload_errors"])
        aggregate["before_bch_payload_bits"] += 640


def export_r8_vector(
    output_dir: Path,
    context: FrameContext,
    result: dict[str, Any],
    artifacts: dict[str, Any],
) -> dict[str, Any]:
    snr_name = f"{context.ebn0_db:.2f}".replace("-", "m").replace(".", "p")
    vector_name = f"ebn0_{snr_name}_seed_{context.seed}_frame_{context.frame_index:06d}"
    vector_dir = output_dir / "hardware_vectors" / vector_name
    vector_dir.mkdir(parents=True, exist_ok=True)

    payload_bytes = bits_to_bytes(context.payload)
    codeword_bytes = bits_to_bytes(artifacts["ldpc_codeword"])
    llr_bytes = b"".join(struct.pack("<i", int(value)) for value in artifacts["llr"])
    (vector_dir / "payload.bin").write_bytes(payload_bytes)
    (vector_dir / "codeword.bin").write_bytes(codeword_bytes)
    (vector_dir / "llr32.bin").write_bytes(llr_bytes)

    payload_out = artifacts.get("payload_out")
    metadata = {
        "name": vector_name,
        "ebn0_db": context.ebn0_db,
        "seed": context.seed,
        "frame_index": context.frame_index,
        "llr_scale": result.get("llr_scale"),
        "payload_bytes": len(payload_bytes),
        "codeword_bytes": len(codeword_bytes),
        "llr_count": len(artifacts["llr"]),
        "llr_storage": "1120 signed values, little-endian int32",
        "payload_sha256": hashlib.sha256(payload_bytes).hexdigest(),
        "codeword_sha256": hashlib.sha256(codeword_bytes).hexdigest(),
        "llr32_sha256": hashlib.sha256(llr_bytes).hexdigest(),
        "ldpc_ok": bool(result["ldpc_ok"]),
        "bch_ran": bool(result["bch_ran"]),
        "bch_ok": None if result["bch_ok"] == "" else bool(result["bch_ok"]),
        "decode_reject": bool(result["decode_reject"]),
        "expected_payload_hex": bits_to_bytes(payload_out).hex() if payload_out is not None else None,
        "expected_result": "success" if result["output_valid"] else "decode_reject",
    }
    (vector_dir / "expected.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return metadata


def self_test(ldpc840: Any, ldpc640: Any, max_iter: int) -> None:
    bch = BCHShortened(n=840, k=640, t=20)
    decoder840 = R8IntegerLayeredDecoder(ldpc840, max_iter=max_iter)
    decoder640 = R8IntegerLayeredDecoder(ldpc640, max_iter=max_iter)
    rng = random.Random(20260911)
    payload = rand_bits(rng, 640)

    bch_codeword = bch.encode_bits(payload)
    ldpc840_codeword = ldpc840.encode_bits(bch_codeword)
    llr840 = [20 if bit == 0 else -20 for bit in ldpc840_codeword]
    decoded840, _hard840, ok840, _it840, _lq840 = decoder840.decode(llr840)
    if not ok840 or decoded840 != bch_codeword:
        raise RuntimeError("Self-test failed: clean R8 LDPC decode mismatch")
    bch_decoded, bch_ok, _corrected = bch.decode_bits(decoded840)
    if not bch_ok or bch_decoded != payload:
        raise RuntimeError("Self-test failed: clean BCH decode mismatch")

    ldpc640_codeword = ldpc640.encode_bits(payload)
    llr640 = [20 if bit == 0 else -20 for bit in ldpc640_codeword]
    decoded640, _hard640, ok640, _it640, _lq640 = decoder640.decode(llr640)
    if not ok640 or decoded640 != payload:
        raise RuntimeError("Self-test failed: clean LDPC(1120,640) decode mismatch")

    error_codeword = list(bch_codeword)
    for position in (3, 61, 117, 239, 511):
        error_codeword[position] ^= 1
    corrected, corrected_ok, corrected_count = bch.decode_bits(error_codeword)
    if not corrected_ok or corrected != payload or corrected_count != 5:
        raise RuntimeError("Self-test failed: BCH five-error correction mismatch")

    print("PASS: R8 software comparison self-test")


def write_plots(summary_path: Path, output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"matplotlib unavailable; plots skipped: {exc}")
        return

    with summary_path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[row["mode"]].append(row)
    for values in grouped.values():
        values.sort(key=lambda item: float(item["ebn0_db"]))

    def plot_metric(column: str, ylabel: str, filename: str) -> None:
        plt.figure(figsize=(10, 6))
        for mode, values in grouped.items():
            x_values = [float(row["ebn0_db"]) for row in values]
            y_values = [max(float(row[column]), 1e-8) for row in values]
            plt.semilogy(
                x_values,
                y_values,
                marker="o",
                linewidth=2.2,
                label=MODE_LABELS.get(mode, mode),
            )
        plt.xlabel("Eb/N0 (dB)")
        plt.ylabel(ylabel)
        plt.grid(True, which="both", linestyle="--", alpha=0.45)
        plt.legend()
        plt.tight_layout()
        path = output_dir / filename
        plt.savefig(path, dpi=180)
        plt.close()
        print(f"Saved {path}")

    plot_metric("effective_output_ber", "Effective output BER", "software_comparison_ber.png")
    plot_metric("fer", "FER", "software_comparison_fer.png")

    bch_mode_for_plot = "bch_ldpc_fallback" if "bch_ldpc_fallback" in grouped else "bch_ldpc_r8"
    bch_rows = grouped.get(bch_mode_for_plot, [])
    if bch_rows:
        x_values = [float(row["ebn0_db"]) for row in bch_rows]
        cleaned = [int(float(row["bch_cleaned_frames"])) for row in bch_rows]
        residual_small = [int(float(row["residual_1_to_20_frames"])) for row in bch_rows]
        residual_large = [int(float(row["residual_gt_20_frames"])) for row in bch_rows]
        width = 0.055 if len(x_values) > 1 else 0.12
        plt.figure(figsize=(10, 6))
        plt.bar([x - width for x in x_values], residual_small, width=width, label="LDPC residual 1..20")
        plt.bar(x_values, residual_large, width=width, label="LDPC residual >20")
        plt.bar([x + width for x in x_values], cleaned, width=width, label="Frames cleaned by BCH")
        plt.xlabel("Eb/N0 (dB)")
        plt.ylabel("Frame count")
        plt.title(f"BCH contribution: {MODE_LABELS.get(bch_mode_for_plot, bch_mode_for_plot)}")
        plt.grid(True, axis="y", linestyle="--", alpha=0.4)
        plt.legend()
        plt.tight_layout()
        path = output_dir / "software_bch_contribution.png"
        plt.savefig(path, dpi=180)
        plt.close()
        print(f"Saved {path}")


def run(args: argparse.Namespace) -> None:
    output_dir = args.out_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "summary.csv"
    frame_path = output_dir / "frames.csv"
    metadata_path = output_dir / "run_metadata.json"

    ldpc840 = make_ldpc_1120_840(args.max_iter)
    ldpc640, ldpc640_name, ldpc640_description = make_ldpc_1120_640(
        args.ldpc640_config.resolve(), args.max_iter
    )
    decoder840 = R8IntegerLayeredDecoder(ldpc840, max_iter=args.max_iter)
    decoder640 = R8IntegerLayeredDecoder(ldpc640, max_iter=args.max_iter)
    bch = BCHShortened(n=840, k=640, t=20)

    if args.self_test:
        self_test(ldpc840, ldpc640, args.max_iter)

    if args.modes:
        modes = [item.strip() for item in args.modes.split(",") if item.strip()]
        invalid_modes = [mode for mode in modes if mode not in AVAILABLE_MODES]
        if invalid_modes:
            raise SystemExit(
                "ERROR: unknown --modes value(s): " + ", ".join(invalid_modes)
                + ". Allowed values: " + ", ".join(AVAILABLE_MODES)
            )
        if len(set(modes)) != len(modes):
            raise SystemExit("ERROR: --modes must not contain duplicates")
    elif args.study == "theory":
        modes = list(THEORY_MODES)
        if args.include_ablation:
            modes.insert(2, "ldpc_640pad")
    elif args.study == "r8":
        modes = list(R8_MODES)
        if args.include_ablation:
            modes.insert(2, "ldpc_640pad")
    else:
        modes = ["uncoded", "ldpc_1120_640", "bch_ldpc_fallback", "bch_ldpc_r8"]
        if args.include_ablation:
            modes.insert(2, "ldpc_640pad")
    bch_modes = [mode for mode in modes if mode.startswith("bch_ldpc_")]
    if args.export_vectors_per_snr > 0 and "bch_ldpc_r8" not in bch_modes:
        raise SystemExit(
            "ERROR: --export-vectors-per-snr requires --study r8 or --study both "
            "because exported vectors carry the expected R8 gate result."
        )
    snr_values = parse_snr_list(args.snr)
    seeds = parse_int_list(args.seeds)

    print(f"Study mode       : {args.study}")
    print(f"Output directory : {output_dir}")
    print(f"Eb/N0 points     : {snr_values}")
    print(f"Seeds            : {seeds}")
    print(f"Frames/seed/point: {args.frames_per_seed}")
    print(f"LLR scale        : {args.llr_scale}")
    print(f"LDPC iterations  : {args.max_iter}")
    print(f"LDPC840          : {ldpc840.summary()}")
    print(f"LDPC640          : {ldpc640.summary()} ({ldpc640_name})")

    metadata = {
        "study": "R8 paired software coding comparison",
        "study_mode": args.study,
        "modes": modes,
        "ebn0_db": snr_values,
        "seeds": seeds,
        "frames_per_seed_per_point": args.frames_per_seed,
        "payload_bits": 640,
        "coded_bits": 1120,
        "coded_rate": 640 / 1120,
        "llr_quantization": "nearest-even then signed-six-bit clip [-32,31]",
        "llr_scale": args.llr_scale,
        "ldpc_alpha": 0.75,
        "ldpc_max_iter": args.max_iter,
        "r8_ldpc": ldpc840.summary(),
        "experimental_ldpc640": ldpc640.summary(),
        "experimental_ldpc640_name": ldpc640_name,
        "experimental_ldpc640_description": ldpc640_description,
        "bch_ldpc_fallback_policy": "BCH runs on every LDPC hard output",
        "bch_ldpc_r8_policy": "BCH runs only when LDPC parity check succeeds",
        "rejected_frame_effective_errors": 640,
        "noise_pairing": "same payload for all modes; same 1120 Gaussian samples for coded modes",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    frame_fields = list(base_result(FrameContext(0.0, 0, 0, [], []), "x").keys()) + ["llr_scale"]
    aggregates: dict[tuple[str, float], dict[str, float]] = defaultdict(make_aggregator)
    paired_outcomes: dict[tuple[float, int, int], dict[str, int]] = defaultdict(dict)
    exported_per_snr: dict[float, int] = defaultdict(int)
    vector_manifest: list[dict[str, Any]] = []

    total_work = len(snr_values) * len(seeds) * args.frames_per_seed
    completed = 0
    with frame_path.open("w", newline="", encoding="utf-8") as frame_handle:
        writer = csv.DictWriter(frame_handle, fieldnames=frame_fields)
        writer.writeheader()

        for ebn0_db in snr_values:
            for seed in seeds:
                rng = make_rng(seed, ebn0_db)
                for frame_index in range(args.frames_per_seed):
                    context = FrameContext(
                        ebn0_db=ebn0_db,
                        seed=seed,
                        frame_index=frame_index,
                        payload=rand_bits(rng, 640),
                        standard_noise=[rng.gauss(0.0, 1.0) for _ in range(1120)],
                    )

                    uncoded = simulate_uncoded(context) if "uncoded" in modes else None
                    ldpc_only = None
                    if "ldpc_1120_640" in modes:
                        ldpc_only, _ldpc_artifacts = simulate_ldpc_only(
                            context,
                            ldpc640,
                            decoder640,
                            args.llr_scale,
                            "ldpc_1120_640",
                            padded=False,
                        )
                    bch_results: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
                    if "bch_ldpc_fallback" in bch_modes:
                        bch_results["bch_ldpc_fallback"] = simulate_bch_ldpc_fallback(
                            context, ldpc840, decoder840, bch, args.llr_scale
                        )
                    if "bch_ldpc_r8" in bch_modes:
                        bch_results["bch_ldpc_r8"] = simulate_bch_ldpc_r8(
                            context, ldpc840, decoder840, bch, args.llr_scale
                        )

                    rows = []
                    if uncoded is not None:
                        rows.append(uncoded)
                    if ldpc_only is not None:
                        rows.append(ldpc_only)
                    if "ldpc_640pad" in modes:
                        ablation, _ablation_artifacts = simulate_ldpc_only(
                            context,
                            ldpc840,
                            decoder840,
                            args.llr_scale,
                            "ldpc_640pad",
                            padded=True,
                        )
                        rows.append(ablation)
                    rows.extend(result for result, _artifacts in bch_results.values())

                    for row in rows:
                        row["llr_scale"] = args.llr_scale
                        writer.writerow(row)
                        update_aggregator(aggregates[(str(row["mode"]), ebn0_db)], row)

                    pair_key = (ebn0_db, seed, frame_index)
                    if ldpc_only is not None:
                        paired_outcomes[pair_key]["ldpc_1120_640"] = int(ldpc_only["frame_error"])
                    for bch_mode, (bch_result, _bch_artifacts) in bch_results.items():
                        paired_outcomes[pair_key][bch_mode] = int(bch_result["frame_error"])

                    if (
                        args.export_vectors_per_snr > 0
                        and exported_per_snr[ebn0_db] < args.export_vectors_per_snr
                    ):
                        bch_result, bch_artifacts = bch_results["bch_ldpc_r8"]
                        vector_manifest.append(
                            export_r8_vector(output_dir, context, bch_result, bch_artifacts)
                        )
                        exported_per_snr[ebn0_db] += 1

                    completed += 1
                    if args.progress_every > 0 and completed % args.progress_every == 0:
                        print(f"Progress: {completed}/{total_work} paired frames")

    summary_fields = [
        "mode",
        "ebn0_db",
        "frames",
        "valid_frames",
        "payload_bits",
        "transmitted_bits",
        "channel_bit_errors",
        "channel_ber",
        "payload_bit_errors_valid",
        "valid_payload_ber",
        "effective_output_bit_errors",
        "effective_output_ber",
        "frame_errors",
        "fer",
        "fer_ci95_low",
        "fer_ci95_high",
        "decode_rejects",
        "accepted_wrong_frames",
        "avg_ldpc_iterations",
        "residual_zero_frames",
        "residual_1_to_20_frames",
        "residual_gt_20_frames",
        "before_bch_payload_ber",
        "after_bch_payload_ber",
        "bch_ran_frames",
        "bch_corrected_errors",
        "bch_cleaned_frames",
    ]
    summary_rows: list[dict[str, Any]] = []
    for mode in modes:
        for ebn0_db in snr_values:
            aggregate = aggregates[(mode, ebn0_db)]
            frames = int(aggregate["frames"])
            frame_errors = int(aggregate["frame_errors"])
            ci_low, ci_high = wilson_interval(frame_errors, frames)
            valid_frames = int(aggregate["valid_frames"])
            payload_bits = int(aggregate["payload_bits"])
            transmitted_bits = int(aggregate["transmitted_bits"])
            before_bits = int(aggregate["before_bch_payload_bits"])
            row = {
                "mode": mode,
                "ebn0_db": ebn0_db,
                "frames": frames,
                "valid_frames": valid_frames,
                "payload_bits": payload_bits,
                "transmitted_bits": transmitted_bits,
                "channel_bit_errors": int(aggregate["channel_bit_errors"]),
                "channel_ber": aggregate["channel_bit_errors"] / transmitted_bits if transmitted_bits else 0.0,
                "payload_bit_errors_valid": int(aggregate["payload_bit_errors_valid"]),
                "valid_payload_ber": aggregate["payload_bit_errors_valid"] / (valid_frames * 640) if valid_frames else 0.0,
                "effective_output_bit_errors": int(aggregate["effective_output_bit_errors"]),
                "effective_output_ber": aggregate["effective_output_bit_errors"] / payload_bits if payload_bits else 0.0,
                "frame_errors": frame_errors,
                "fer": frame_errors / frames if frames else 0.0,
                "fer_ci95_low": ci_low,
                "fer_ci95_high": ci_high,
                "decode_rejects": int(aggregate["decode_rejects"]),
                "accepted_wrong_frames": int(aggregate["accepted_wrong_frames"]),
                "avg_ldpc_iterations": aggregate["ldpc_iterations"] / frames if frames else 0.0,
                "residual_zero_frames": int(aggregate["residual_zero_frames"]),
                "residual_1_to_20_frames": int(aggregate["residual_1_to_20_frames"]),
                "residual_gt_20_frames": int(aggregate["residual_gt_20_frames"]),
                "before_bch_payload_ber": aggregate["before_bch_payload_errors"] / before_bits if before_bits else 0.0,
                "after_bch_payload_ber": aggregate["after_bch_payload_errors"] / aggregate["after_bch_payload_bits"] if aggregate["after_bch_payload_bits"] else 0.0,
                "bch_ran_frames": int(aggregate["bch_ran_frames"]),
                "bch_corrected_errors": int(aggregate["bch_corrected_errors"]),
                "bch_cleaned_frames": int(aggregate["bch_cleaned_frames"]),
            }
            summary_rows.append(row)
            print(
                f"{mode:15s} Eb/N0={ebn0_db:4.1f} frames={frames:6d} "
                f"effBER={row['effective_output_ber']:.3e} FER={row['fer']:.3e} "
                f"reject={row['decode_rejects']} wrongAccepted={row['accepted_wrong_frames']}"
            )

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)

    paired_fields = [
        "ebn0_db",
        "frames",
        "both_pass",
        "both_fail",
        "ldpc_fail_bch_pass",
        "ldpc_pass_bch_fail",
    ]
    for bch_mode in bch_modes:
        if "ldpc_1120_640" not in modes:
            continue
        paired_rows: list[dict[str, Any]] = []
        for ebn0_db in snr_values:
            counters = defaultdict(int)
            for (row_snr, _seed, _frame), outcome in paired_outcomes.items():
                if row_snr != ebn0_db:
                    continue
                if "ldpc_1120_640" not in outcome or bch_mode not in outcome:
                    continue
                ldpc_fail = outcome["ldpc_1120_640"]
                bch_fail = outcome[bch_mode]
                if not ldpc_fail and not bch_fail:
                    counters["both_pass"] += 1
                elif ldpc_fail and bch_fail:
                    counters["both_fail"] += 1
                elif ldpc_fail and not bch_fail:
                    counters["ldpc_fail_bch_pass"] += 1
                else:
                    counters["ldpc_pass_bch_fail"] += 1
            paired_rows.append(
                {
                    "ebn0_db": ebn0_db,
                    "frames": sum(counters.values()),
                    **{field: counters[field] for field in paired_fields[2:]},
                }
            )
        paired_path = output_dir / f"paired_ldpc_vs_{bch_mode}.csv"
        with paired_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=paired_fields)
            writer.writeheader()
            writer.writerows(paired_rows)

    if vector_manifest:
        (output_dir / "hardware_vectors" / "manifest.json").write_text(
            json.dumps(vector_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )

    write_plots(summary_path, output_dir)
    print(f"PASS: software comparison completed: {summary_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="R8 paired coding comparison")
    parser.add_argument("--snr", default="3.5,4.0,4.2,4.4,4.6,4.8,5.0")
    parser.add_argument("--seeds", default="101,202,303,404,505")
    parser.add_argument("--frames-per-seed", type=int, default=200)
    parser.add_argument("--llr-scale", type=float, default=1.0)
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--ldpc640-config", type=Path, default=DEFAULT_LDPC640_CONFIG)
    parser.add_argument(
        "--study",
        choices=("theory", "r8", "both"),
        default="theory",
        help="theory: BCH fallback; r8: deployed RTL gate; both: run both BCH policies",
    )
    parser.add_argument("--include-ablation", action="store_true")
    parser.add_argument(
        "--modes",
        default="",
        help=(
            "Comma-separated explicit mode subset. Overrides --study and "
            "--include-ablation; allowed: " + ", ".join(AVAILABLE_MODES)
        ),
    )
    parser.add_argument("--export-vectors-per-snr", type=int, default=0)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.frames_per_seed <= 0:
        parser.error("--frames-per-seed must be positive")
    if args.llr_scale <= 0.0 or not math.isfinite(args.llr_scale):
        parser.error("--llr-scale must be finite and positive")
    if args.max_iter <= 0:
        parser.error("--max-iter must be positive")
    if args.export_vectors_per_snr < 0:
        parser.error("--export-vectors-per-snr cannot be negative")
    run(args)


if __name__ == "__main__":
    main()
