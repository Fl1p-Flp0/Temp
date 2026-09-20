[CmdletBinding()]
param(
    [ValidateRange(1, 100000000)]
    [int]$Frames = 5000,

    [int[]]$Seeds = @(20260920, 20260921, 20260922),

    [string]$Snr = "2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2"
)

$ErrorActionPreference = "Stop"
$study = Join-Path $PSScriptRoot "model\python\Secure_ECC\link_sim\bch_ldpc_iterative_study.py"
if (-not (Test-Path -LiteralPath $study)) {
    throw "Simulation script not found: $study"
}

if (Get-Command py -ErrorAction SilentlyContinue) {
    $python = "py"
    $pythonPrefix = @("-3")
} elseif (Get-Command python -ErrorAction SilentlyContinue) {
    $python = "python"
    $pythonPrefix = @()
} else {
    throw "Python 3 was not found. Install Python 3.9 or newer, then rerun."
}

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$root = Join-Path $PSScriptRoot "results\r9_768b_onepass_$stamp"
foreach ($seed in $Seeds) {
    $outDir = Join-Path $root "seed_$seed"
    $arguments = @(
        $pythonPrefix
        $study
        "--snr", $Snr
        "--frames", $Frames
        "--seed", $seed
        "--max-iter", "10"
        "--application-bits", "768"
        "--modes", "ldpc_only,onepass"
        "--out-dir", $outDir
    )

    Write-Host "Running paired LDPC-only versus BCH one-pass, seed=$seed -> $outDir"
    & $python @arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Simulation failed for seed=$seed with exit code $LASTEXITCODE"
    }
}

Write-Host "PASS: paired LDPC-only versus BCH one-pass study completed"
Write-Host "RESULTS=$root"
