# CPU Benchmark: AES-BCH-LDPC Processing Chain

This repository contains a self-contained host-CPU benchmark for a 768-bit
useful-payload processing chain:

```text
TX: 96-byte plaintext -> AES-128 x6 -> BCH(840,800,t=4) -> LDPC(1120,840)
RX: signed 6-bit LLR -> LDPC(1120,840) -> BCH(840,800,t=4) -> AES-128 x6
```

The BCH and LDPC masks are embedded in `cpu_benchmark/rtl_masks.py`. No FPGA
project, synthesis tool, driver, DMA setup, or third-party Python package is
required.

## Measurements

| Measurement | Channel condition | Purpose |
| --- | --- | --- |
| TX | Clean processing path | Noise is injected only after TX has emitted its codeword, so it cannot change TX encoding work. |
| RX clean | Deterministic signed 6-bit LLR values | Measures receive processing capability without channel-induced decoder variation. |
| RX AWGN | BPSK AWGN converted to signed 6-bit LLR values | Measures realistic LDPC convergence, BCH behavior, and fail-closed rejects. |

Therefore, TX has one benchmark only, while RX has both clean and AWGN
benchmarks. Each frame carries 768 useful bits.

The primary latency is the **median of the per-repetition average latency per
frame**. Mean and population standard deviation are recorded only to show
measurement repeatability.

```text
payload throughput (Mbit/s) = 768 / latency_us_per_frame
accepted goodput (Mbit/s) = payload throughput * accepted_frames / total_frames
```

For AWGN, use accepted goodput with raw processing throughput because a
fail-closed RX may reject a frame rather than output an incorrect plaintext.

## Quick Start on Another CPU

Requirements:

- Python 3.9 or newer.
- Windows PowerShell 5.1 or PowerShell 7 for the convenience runner.
- No `pip install` step is needed; only the Python standard library is used.

Clone and run the complete standard set:

```powershell
git clone https://github.com/Fl1p-Flp0/Temp.git
cd Temp
Set-ExecutionPolicy -Scope Process Bypass
.\cpu_benchmark\run_cpu_benchmarks.ps1
```

The defaults are intentionally consistent for a general CPU comparison:

```text
TX:       1,000 frames x 50 repetitions, 50 warm-up frames/repetition
RX clean:   100 frames x 50 repetitions,  5 warm-up frames/repetition
RX AWGN:    100 frames x 50 repetitions,  5 warm-up frames/repetition
Eb/N0: 4.0 dB, LLR scale: 1.0, seed: 20260922
```

Run a short smoke test first when useful:

```powershell
.\cpu_benchmark\run_cpu_benchmarks.ps1 `
  -TxFrames 10 -RxFrames 2 -TxWarmup 2 -RxWarmup 1 -Repetitions 1
```

Change the AWGN point or output directory:

```powershell
.\cpu_benchmark\run_cpu_benchmarks.ps1 `
  -EbN0Db 4.0 `
  -LlrScale 1.0 `
  -Seed 20260922 `
  -OutDir build\cpu_benchmark_other_cpu
```

Rerun only one portion when necessary:

```powershell
.\cpu_benchmark\run_cpu_benchmarks.ps1 -Mode tx
.\cpu_benchmark\run_cpu_benchmarks.ps1 -Mode rx-clean
.\cpu_benchmark\run_cpu_benchmarks.ps1 -Mode rx-awgn -EbN0Db 4.0
```

## Output Files

The default root is `build/cpu_benchmark`:

```text
tx_clean/
  tx_python_summary.csv
  tx_python_runs.csv
  tx_python_summary.json
  reference_tx_codeword.hex

rx_clean/
  rx_python_summary.csv
  rx_python_runs.csv
  rx_python_summary.json
  llr6_clean_corpus.bin

rx_awgn_4p0db/
  rx_python_summary.csv
  rx_python_runs.csv
  rx_python_summary.json
  llr6_awgn_corpus.bin
```

For result aggregation, keep the three summary CSV and JSON files, all three
`*_runs.csv` files, and the AWGN corpus. The corpus contains one signed 6-bit
LLR value per coded bit and can be replayed by compatible software, RV32, or
FPGA receive benchmarks.

## Direct Python Commands

The benchmark also works on Linux and macOS without PowerShell:

```text
python3 cpu_benchmark/benchmark_tx_python.py \
  --frames 1000 --warmup 50 --repetitions 50 \
  --out-dir build/cpu_benchmark/tx_clean

python3 cpu_benchmark/benchmark_rx_python.py \
  --scenario clean --frames 100 --warmup 5 --repetitions 50 \
  --out-dir build/cpu_benchmark/rx_clean

python3 cpu_benchmark/benchmark_rx_python.py \
  --scenario awgn --ebn0-db 4.0 --llr-scale 1.0 --seed 20260922 \
  --frames 100 --warmup 5 --repetitions 50 \
  --write-corpus build/cpu_benchmark/rx_awgn_4p0db/llr6_awgn_corpus.bin \
  --out-dir build/cpu_benchmark/rx_awgn_4p0db
```

For a fair CPU, RV32, and FPGA comparison, preserve payload size, RX corpus,
decoder configuration, Eb/N0, LLR scale, and primary median statistic. Do not
compare CPU core-only latency directly with board transaction time that also
contains DMA, Linux scheduling, and driver overhead.
