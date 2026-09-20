# R9 BCH-Aided LDPC Simulation

This repository is a portable software experiment for the R9 candidate:

```text
768 application bits + 32 fixed zero bits
-> BCH(840,800,t=4)
-> QC-LDPC(1120,840)
-> AWGN with signed 6-bit LLR
```

The runner compares four modes for each frame using the same application
payload and the same standard-normal AWGN samples:

1. `LDPC(1120,840) only`: 768 application bits + 72 zero pad bits -> LDPC.
2. `LDPC terminal (BCH bypass)`: BCH is present at the transmitter but not at
   the receiver. This is a diagnostic baseline, not an LDPC-only system.
3. `LDPC -> BCH one-pass`: ten LDPC iterations then one BCH decode.
4. `BCH-aided LDPC (5+5)`: five LDPC iterations -> BCH decode -> full 840-bit
   BCH decision feedback -> five fresh LDPC iterations -> final BCH decode.

The LDPC-only and BCH-LDPC systems both carry 768 useful bits in a 1120-bit
frame, giving the same useful payload rate of `768 / 1120 = 0.685714`.

## Requirements

- Windows PowerShell 5.1 or PowerShell 7.
- Python 3.9 or newer. The `py -3` launcher is preferred; `python` is also
  supported by the runner.
- `matplotlib` for PNG plots.

Install the Python dependency once:

```powershell
py -3 -m pip install -r requirements.txt
```

## Gamma Pilot

Run the default pilot sweep. It evaluates gamma values 2, 4, 6 and 8 over
eight Eb/N0 points with 500 frames per point:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\run_r9_768b_gamma.ps1
```

To run a longer experiment after selecting a candidate gamma:

```powershell
.\run_r9_768b_gamma.ps1 -Frames 5000 -Gammas 4
```

## Outputs

Each gamma directory in `results/` contains:

```text
summary.csv
run_metadata.json
bch_ldpc_t4_effective_ber.png
bch_ldpc_t4_fer.png
```

Send the complete `results/r9_768b_gamma_<timestamp>/` directory back for
aggregation. The CSV contains raw-decision BER, fail-closed effective-output
BER, FER, BCH success and rejection counts, and BCH feedback statistics.

## Reproducibility Notes

- BCH feedback uses all 840 corrected BCH decisions, including BCH parity
  positions 800 through 839.
- Feedback has a configurable signed-six-bit magnitude named `gamma`.
- The second LDPC pass restarts its check-node messages. This is a
  hardware-friendly two-pass candidate, not a bit-exact reproduction of a
  published extrinsic-information decoder.
