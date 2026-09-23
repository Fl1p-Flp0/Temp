[CmdletBinding()]
param(
    [ValidateSet("all", "tx", "rx-clean", "rx-awgn")]
    [string]$Mode = "all",
    [ValidateRange(1, 100000000)]
    [int]$TxFrames = 1000,
    [ValidateRange(1, 100000000)]
    [int]$RxFrames = 100,
    [ValidateRange(0, 100000000)]
    [int]$TxWarmup = 50,
    [ValidateRange(0, 100000000)]
    [int]$RxWarmup = 5,
    [ValidateRange(1, 100000)]
    [int]$Repetitions = 50,
    [double]$EbN0Db = 4.0,
    [ValidateRange(0.000001, 1000.0)]
    [double]$LlrScale = 1.0,
    [long]$Seed = 20260922,
    [string]$OutDir = "build\\cpu_benchmark"
)

$ErrorActionPreference = "Stop"
if (Get-Command py -ErrorAction SilentlyContinue) {
    $python = "py"
    $pythonPrefix = @("-3")
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $python = "python"
    $pythonPrefix = @()
} else {
    throw "Python 3 was not found. Install Python 3.9 or newer, then rerun."
}

$txScript = Join-Path $PSScriptRoot "benchmark_tx_python.py"
$rxScript = Join-Path $PSScriptRoot "benchmark_rx_python.py"
foreach ($path in @($txScript, $rxScript)) {
    if (-not (Test-Path -LiteralPath $path)) { throw "Benchmark script not found: $path" }
}

$outputRoot = [System.IO.Path]::GetFullPath((Join-Path (Get-Location) $OutDir))
New-Item -ItemType Directory -Force -Path $outputRoot | Out-Null

function Run-TxClean {
    $output = Join-Path $outputRoot "tx_clean"
    Write-Host "[1/3] TX clean processing path -> $output"
    & $python @pythonPrefix $txScript `
        "--frames" $TxFrames "--warmup" $TxWarmup "--repetitions" $Repetitions `
        "--out-dir" $output "--emit-codeword" (Join-Path $output "reference_tx_codeword.hex")
    if ($LASTEXITCODE -ne 0) { throw "TX benchmark failed with exit code $LASTEXITCODE" }
}

function Run-RxClean {
    $output = Join-Path $outputRoot "rx_clean"
    Write-Host "[2/3] RX clean signed-6-bit LLR path -> $output"
    & $python @pythonPrefix $rxScript `
        "--scenario" "clean" "--frames" $RxFrames "--warmup" $RxWarmup "--repetitions" $Repetitions `
        "--out-dir" $output "--write-corpus" (Join-Path $output "llr6_clean_corpus.bin")
    if ($LASTEXITCODE -ne 0) { throw "Clean RX benchmark failed with exit code $LASTEXITCODE" }
}

function Run-RxAwgn {
    $snrTag = $EbN0Db.ToString("0.0", [System.Globalization.CultureInfo]::InvariantCulture).Replace(".", "p")
    $output = Join-Path $outputRoot "rx_awgn_$($snrTag)db"
    Write-Host "[3/3] RX AWGN path, Eb/N0=$EbN0Db dB -> $output"
    & $python @pythonPrefix $rxScript `
        "--scenario" "awgn" `
        "--ebn0-db" $EbN0Db.ToString([System.Globalization.CultureInfo]::InvariantCulture) `
        "--llr-scale" $LlrScale.ToString([System.Globalization.CultureInfo]::InvariantCulture) `
        "--seed" $Seed "--frames" $RxFrames "--warmup" $RxWarmup "--repetitions" $Repetitions `
        "--out-dir" $output "--write-corpus" (Join-Path $output "llr6_awgn_corpus.bin")
    if ($LASTEXITCODE -ne 0) { throw "AWGN RX benchmark failed with exit code $LASTEXITCODE" }
}

switch ($Mode) {
    "all" { Run-TxClean; Run-RxClean; Run-RxAwgn }
    "tx" { Run-TxClean }
    "rx-clean" { Run-RxClean }
    "rx-awgn" { Run-RxAwgn }
}

Write-Host "PASS: requested CPU benchmark(s) completed"
Write-Host "RESULTS=$outputRoot"
