#!/usr/bin/env python3
"""Paired BCH-aided two-pass LDPC experiment for the R9 candidate.

Every requested receiver mode uses the same application payload and AWGN
sample for a given frame. The study includes a fair system baseline:

    application payload + zero padding -> LDPC(1120,840)

and the BCH-aided candidate:

    application payload + zero padding -> BCH(840,800,t=4) -> LDPC(1120,840)

Receiver modes:

* terminal: ten LDPC iterations and BCH bypassed;
* onepass: ten LDPC iterations followed by one BCH decode;
* iterative: five LDPC iterations, BCH decode, full 840-bit BCH feedback,
  five fresh LDPC iterations, then a final BCH decode.

The second LDPC pass deliberately restarts check-node messages. This is a
hardware-friendly two-pass candidate, not a bit-exact reproduction of a
published extrinsic-information decoder.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from r8_software_coding_comparison import R8IntegerLayeredDecoder, channel_llr6
from secure_ecc_link_sim import BCHShortened, count_errors, import_ldpc_model, parse_snr_list, rand_bits


MODE_LABELS = {
    "ldpc_only": "LDPC(1120,840) only",
    "terminal": "LDPC terminal (BCH bypass)",
    "onepass": "LDPC -> BCH one-pass",
    "iterative": "BCH-aided LDPC (5+5)",
}

PLOT_STYLES = {
    MODE_LABELS["ldpc_only"]: ("#222222", "*"),
    MODE_LABELS["terminal"]: ("#7f8c8d", "x"),
    MODE_LABELS["onepass"]: ("#2471a3", "s"),
    MODE_LABELS["iterative"]: ("#c0392b", "o"),
}


def make_ldpc(max_iter: int) -> Any:
    model = import_ldpc_model()
    cfg = dict(model.PRESETS["1120_840"])
    return model.QCLDPC(cfg["pcm"], int(cfg["z"]), 0.75, max_iter)


def rng_for(seed: int, snr: float) -> random.Random:
    return random.Random(f"bch-iter-r9:{seed}:{snr:.8f}")


def upper95_zero(trials: int) -> float:
    return 1.0 - 0.05 ** (1.0 / trials) if trials else 1.0


@dataclass(frozen=True)
class FrameContext:
    application: list[int]
    bch_message: list[int]
    bch_codeword: list[int]
    standard_noise: list[float]
    llr: list[int]
    channel_errors: int
    ebn0_db: float


def make_frame_context(
    snr: float,
    rng: random.Random,
    code: Any,
    bch: BCHShortened,
    application_bits: int,
    scale: float,
) -> FrameContext:
    application = rand_bits(rng, application_bits)
    bch_message = application + [0] * (bch.k - application_bits)
    bch_codeword = bch.encode_bits(bch_message)
    ldpc_codeword = code.encode_bits(bch_codeword)
    noise = [rng.gauss(0.0, 1.0) for _ in range(1120)]
    llr, channel_errors = channel_llr6(
        ldpc_codeword, snr, application_bits / 1120, noise, scale
    )
    return FrameContext(application, bch_message, bch_codeword, noise, llr, channel_errors, snr)


def apply_bch_feedback(
    llr: list[int],
    first_ldpc_word: list[int],
    corrected_bch_word: list[int],
    gamma: int,
    feedback_mode: str,
) -> tuple[list[int], int]:
    """Add finite BCH evidence to all 840 systematic LDPC positions.

    ``all`` sends signed a-priori evidence for the complete corrected BCH
    word. ``changed`` updates only positions altered by BCH. Both modes retain
    all 840 BCH positions, including the BCH parity bits at indices 800..839.
    """
    if len(first_ldpc_word) != 840 or len(corrected_bch_word) != 840:
        raise ValueError("BCH feedback requires a complete 840-bit codeword")
    feedback = list(llr)
    applied = 0
    for index, corrected_bit in enumerate(corrected_bch_word):
        if feedback_mode == "changed" and first_ldpc_word[index] == corrected_bit:
            continue
        evidence = gamma if corrected_bit == 0 else -gamma
        feedback[index] = max(-32, min(31, feedback[index] + evidence))
        applied += 1
    return feedback, applied


def finalize_bch_result(
    context: FrameContext,
    output_message: list[int],
    bch_ok: bool,
    stats: dict[str, int],
) -> dict[str, int]:
    application_bits = len(context.application)
    raw_errors = count_errors(context.application, output_message[:application_bits])
    padding_ok = all(bit == 0 for bit in output_message[application_bits:])
    accepted = bool(bch_ok and padding_ok)
    stats.update(
        bits=application_bits,
        raw_decision_errors=raw_errors,
        accepted_bit_errors=raw_errors if accepted else 0,
        effective_output_errors=raw_errors if accepted else application_bits,
        frame_errors=int((not accepted) or raw_errors > 0),
        rejects=int(not accepted),
        accepted_wrong_frames=int(accepted and raw_errors > 0),
        padding_failures=int(not padding_ok),
    )
    return stats


def run_ldpc_only(
    context: FrameContext,
    code: Any,
    decoder: R8IntegerLayeredDecoder,
    scale: float,
) -> dict[str, int]:
    """LDPC-only system with the same useful payload and 1120-bit frame size.

    For example, a 768-bit AES payload is followed by 72 fixed zero bits to
    form the 840-bit LDPC message. The BCH-aided path instead uses 32 fixed
    zeros followed by 40 BCH parity bits. Both systems therefore carry the
    same useful payload in the same 1120 transmitted bits.
    """
    application_bits = len(context.application)
    message = context.application + [0] * (840 - application_bits)
    codeword = code.encode_bits(message)
    llr, channel_errors = channel_llr6(
        codeword,
        context.ebn0_db,
        application_bits / 1120,
        context.standard_noise,
        scale,
    )
    output, _hard, syndrome_ok, used, _state = decoder.decode(llr)
    raw_errors = count_errors(context.application, output[:application_bits])
    padding_ok = all(bit == 0 for bit in output[application_bits:840])
    accepted = bool(syndrome_ok and padding_ok)
    return {
        "bits": application_bits,
        "raw_decision_errors": raw_errors,
        "accepted_bit_errors": raw_errors if accepted else 0,
        "effective_output_errors": raw_errors if accepted else application_bits,
        "frame_errors": int((not accepted) or raw_errors > 0),
        "rejects": int(not accepted),
        "accepted_wrong_frames": int(accepted and raw_errors > 0),
        "padding_failures": int(not padding_ok),
        "raw_errors": channel_errors,
        "iterations": used,
        "syndrome_fail": int(not syndrome_ok),
        "feedback_frames": 0,
        "feedback_bits": 0,
        "bch_attempts": 0,
        "bch_successes": 0,
    }


def run_terminal(context: FrameContext, decoder: R8IntegerLayeredDecoder) -> dict[str, int]:
    output, _hard, syndrome_ok, used, _state = decoder.decode(context.llr)
    application_bits = len(context.application)
    raw_errors = count_errors(context.application, output[:application_bits])
    return {
        "bits": application_bits,
        "raw_decision_errors": raw_errors,
        "accepted_bit_errors": raw_errors,
        "effective_output_errors": raw_errors,
        "frame_errors": int(bool(raw_errors)),
        "rejects": 0,
        "accepted_wrong_frames": int(bool(raw_errors)),
        "padding_failures": 0,
        "raw_errors": context.channel_errors,
        "iterations": used,
        "syndrome_fail": int(not syndrome_ok),
        "feedback_frames": 0,
        "feedback_bits": 0,
        "bch_attempts": 0,
        "bch_successes": 0,
    }


def run_onepass(
    context: FrameContext,
    decoder: R8IntegerLayeredDecoder,
    bch: BCHShortened,
) -> dict[str, int]:
    word, _hard, syndrome_ok, used, _state = decoder.decode(context.llr)
    output, _corrected, bch_ok, _count = bch.decode_codeword(word[:840])
    return finalize_bch_result(context, output, bch_ok, {
        "raw_errors": context.channel_errors,
        "iterations": used,
        "syndrome_fail": int(not syndrome_ok),
        "feedback_frames": 0,
        "feedback_bits": 0,
        "bch_attempts": 1,
        "bch_successes": int(bch_ok),
    })


def run_iterative(
    context: FrameContext,
    decoder1: R8IntegerLayeredDecoder,
    decoder2: R8IntegerLayeredDecoder,
    bch: BCHShortened,
    gamma: int,
    feedback_mode: str,
) -> dict[str, int]:
    first, _hard, first_syndrome_ok, used1, _state = decoder1.decode(context.llr)
    first_word = first[:840]
    _first_output, corrected_word, bch_ok1, _count1 = bch.decode_codeword(first_word)

    feedback_bits = 0
    second_llr = context.llr
    if bch_ok1:
        second_llr, feedback_bits = apply_bch_feedback(
            context.llr, first_word, corrected_word, gamma, feedback_mode
        )

    second, _hard, second_syndrome_ok, used2, _state = decoder2.decode(second_llr)
    output, _corrected, bch_ok2, _count2 = bch.decode_codeword(second[:840])
    return finalize_bch_result(context, output, bch_ok2, {
        "raw_errors": context.channel_errors,
        "iterations": used1 + used2,
        "syndrome_fail": int(not first_syndrome_ok) + int(not second_syndrome_ok),
        "feedback_frames": int(bch_ok1),
        "feedback_bits": feedback_bits,
        "bch_attempts": 2,
        "bch_successes": int(bch_ok1) + int(bch_ok2),
    })


def summarize(name: str, rate: float, snr: float, stats: dict[str, int]) -> dict[str, Any]:
    frames = stats["frames"]
    total_bits = stats["bits"]
    raw_errors = stats["raw_decision_errors"]
    effective_errors = stats["effective_output_errors"]
    frame_errors = stats["frame_errors"]
    return {
        "mode": name,
        "ebn0_db": snr,
        "rate": rate,
        "frames": frames,
        "payload_bits_per_frame": total_bits // frames,
        "total_payload_bits": total_bits,
        "channel_ber": stats["raw_errors"] / (1120 * frames),
        "raw_decision_bit_errors": raw_errors,
        "raw_decision_ber": raw_errors / total_bits,
        "effective_output_bit_errors": effective_errors,
        "effective_output_ber": effective_errors / total_bits,
        "effective_ber_zero_upper95": upper95_zero(total_bits) if not effective_errors else "",
        "frame_errors": frame_errors,
        "fer": frame_errors / frames,
        "fer_zero_upper95": upper95_zero(frames) if not frame_errors else "",
        "avg_ldpc_iterations": stats["iterations"] / frames,
        "ldpc_syndrome_failures": stats["syndrome_fail"],
        "feedback_frames": stats["feedback_frames"],
        "feedback_bits": stats["feedback_bits"],
        "bch_attempts": stats["bch_attempts"],
        "bch_successes": stats["bch_successes"],
        "decode_rejects": stats["rejects"],
        "accepted_wrong_frames": stats["accepted_wrong_frames"],
        "padding_failures": stats["padding_failures"],
    }


def plot(rows: list[dict[str, Any]], out_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"matplotlib unavailable: {exc}")
        return

    def values_for(mode: str) -> list[dict[str, Any]]:
        return sorted(
            (row for row in rows if row["mode"] == mode),
            key=lambda row: float(row["ebn0_db"]),
        )

    fig, ax = plt.subplots(figsize=(10.5, 7))
    for mode, (color, marker) in PLOT_STYLES.items():
        values = values_for(mode)
        if not values:
            continue
        x = [float(row["ebn0_db"]) for row in values]
        observed = [float(row["effective_output_ber"]) for row in values]
        y = [
            value if value else float(row["effective_ber_zero_upper95"])
            for value, row in zip(observed, values)
        ]
        ax.semilogy(x, y, color=color, marker=marker, linewidth=2.2,
                    label=f"{mode}, R={float(values[0]['rate']):.3f}")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("Effective output BER")
    ax.set_title("BCH(840,800,t=4)-aided two-pass LDPC decoding")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "bch_ldpc_t4_effective_ber.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 7))
    for mode, (color, marker) in PLOT_STYLES.items():
        values = values_for(mode)
        if not values:
            continue
        x = [float(row["ebn0_db"]) for row in values]
        observed = [float(row["fer"]) for row in values]
        y = [
            value if value else float(row["fer_zero_upper95"])
            for value, row in zip(observed, values)
        ]
        ax.semilogy(x, y, color=color, marker=marker, linewidth=2.2,
                    label=f"{mode}, R={float(values[0]['rate']):.3f}")
    ax.set_xlabel("Eb/N0 (dB)")
    ax.set_ylabel("FER")
    ax.set_title("BCH(840,800,t=4)-aided two-pass LDPC decoding")
    ax.grid(True, which="both", linestyle=":", alpha=0.6)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "bch_ldpc_t4_fer.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snr", default="3.0,3.2,3.5,3.7,4.0,4.2,4.4,4.6")
    parser.add_argument("--frames", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260920)
    parser.add_argument("--max-iter", type=int, default=10)
    parser.add_argument("--first-iter", type=int, default=5)
    parser.add_argument("--llr-scale", type=float, default=1.0)
    parser.add_argument("--application-bits", type=int, default=640,
                        help="Useful payload bits; the remainder of the 800-bit BCH message is zero padding")
    parser.add_argument("--feedback-gamma", type=int, default=4,
                        help="Signed-six-bit BCH feedback magnitude added to channel LLRs")
    parser.add_argument("--feedback-mode", choices=("all", "changed"), default="all",
                        help="Feedback all 840 BCH decisions or only BCH-corrected positions")
    parser.add_argument("--modes", default="ldpc_only,terminal,onepass,iterative",
                        help="Comma-separated: ldpc_only, terminal, onepass, iterative; aliases ldpc=ldpc_only, t4=iterative")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.frames < 1 or args.max_iter < 2 or args.first_iter < 1 or args.first_iter >= args.max_iter:
        parser.error("Require positive frames and 0 < first-iter < max-iter")
    if args.llr_scale <= 0 or not math.isfinite(args.llr_scale):
        parser.error("llr-scale must be finite and positive")
    if not 1 <= args.application_bits <= 800:
        parser.error("application-bits must be in 1..800")
    if not 1 <= args.feedback_gamma <= 31:
        parser.error("feedback-gamma must be in 1..31")

    aliases = {"ldpc": "ldpc_only", "bch_terminal": "terminal", "t4": "iterative"}
    requested_modes: list[str] = []
    for item in (part.strip().lower() for part in args.modes.split(",")):
        if not item:
            continue
        mode = aliases.get(item, item)
        if mode not in MODE_LABELS:
            parser.error("modes must be ldpc_only, terminal, onepass, iterative; ldpc and t4 are aliases")
        if mode not in requested_modes:
            requested_modes.append(mode)
    if not requested_modes:
        parser.error("At least one mode is required")

    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    code = make_ldpc(args.max_iter)
    full = R8IntegerLayeredDecoder(code, args.max_iter)
    first = R8IntegerLayeredDecoder(code, args.first_iter)
    second = R8IntegerLayeredDecoder(code, args.max_iter - args.first_iter)
    bch = BCHShortened(n=840, k=800, t=4)

    print(f"LDPC: {code.summary()}")
    print(f"BCH: BCH(840,800,t=4), generator degree={bch.g_degree}")
    print(f"Schedule: {args.first_iter}+{args.max_iter - args.first_iter} LDPC iterations")
    print(f"Application payload: {args.application_bits} bits; fixed BCH padding: {800 - args.application_bits} bits")
    print(f"Feedback: {args.feedback_mode}, gamma={args.feedback_gamma}, full 840-bit BCH codeword")
    print(f"Modes: {','.join(requested_modes)}")

    rows: list[dict[str, Any]] = []
    snrs = parse_snr_list(args.snr)
    for snr in snrs:
        stats_by_mode = {mode: defaultdict(int) for mode in requested_modes}
        point_rng = rng_for(args.seed, snr)
        for _ in range(args.frames):
            context = make_frame_context(
                snr, point_rng, code, bch, args.application_bits, args.llr_scale
            )
            runners = {
                "ldpc_only": lambda: run_ldpc_only(context, code, full, args.llr_scale),
                "terminal": lambda: run_terminal(context, full),
                "onepass": lambda: run_onepass(context, full, bch),
                "iterative": lambda: run_iterative(
                    context, first, second, bch, args.feedback_gamma, args.feedback_mode
                ),
            }
            for mode in requested_modes:
                for key, value in runners[mode]().items():
                    stats_by_mode[mode][key] += value
                stats_by_mode[mode]["frames"] += 1

        point_rows: list[dict[str, Any]] = []
        for mode in requested_modes:
            row = summarize(
                MODE_LABELS[mode], args.application_bits / 1120, snr, stats_by_mode[mode]
            )
            rows.append(row)
            point_rows.append(row)
        print(" | ".join(
            f"{row['mode']}: effBER={row['effective_output_ber']:.3e} FER={row['fer']:.3e}"
            for row in point_rows
        ))

    with (out_dir / "summary.csv").open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    (out_dir / "run_metadata.json").write_text(json.dumps({
        "study": "Paired BCH-aided two-pass LDPC candidate",
        "not_paper_bit_exact_reproduction": True,
        "inner_ldpc": "QC-LDPC(1120,840), signed-6-bit NMS",
        "outer_bch": "BCH(840,800,t=4)",
        "ldpc_only_baseline": "application payload + zero padding to 840 bits -> LDPC(1120,840)",
        "schedule": [args.first_iter, args.max_iter - args.first_iter],
        "modes": requested_modes,
        "application_bits": args.application_bits,
        "fixed_zero_padding_bits": 800 - args.application_bits,
        "feedback": {
            "mode": args.feedback_mode,
            "gamma": args.feedback_gamma,
            "codeword_bits": 840,
        },
        "common_payload_and_noise_per_frame": True,
        "frames": args.frames,
        "snr": snrs,
    }, indent=2) + "\n", encoding="utf-8")
    plot(rows, out_dir)
    print(f"PASS: saved {out_dir / 'summary.csv'}")


if __name__ == "__main__":
    main()
