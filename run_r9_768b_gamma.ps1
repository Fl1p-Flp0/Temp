[CmdletBinding()]
param(
    [ValidateRange(1, 100000000)]
    [int]$Frames = 500,

    [int]$Seed = 20260920,

    [string]$Snr = "2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2",

    [int[]]$Gammas = @(2, 4, 6, 8)
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
$root = Join-Path $PSScriptRoot "results\r9_768b_gamma_$stamp"
foreach ($gamma in $Gammas) {
    if ($gamma -lt 1 -or $gamma -gt 31) {
        throw "Each gamma must be in 1..31. Invalid value: $gamma"
    }

    $outDir = Join-Path $root "gamma_$gamma"
    $arguments = @(
        $pythonPrefix
        $study
        "--snr", $Snr
        "--frames", $Frames
        "--seed", $Seed
        "--max-iter", "10"
        "--first-iter", "5"
        "--application-bits", "768"
        "--feedback-gamma", $gamma
        "--feedback-mode", "all"
        "--modes", "ldpc_only,terminal,onepass,iterative"
        "--out-dir", $outDir
    )

    Write-Host "Running gamma=$gamma -> $outDir"
    & $python @arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Simulation failed for gamma=$gamma with exit code $LASTEXITCODE"
    }
}

Write-Host "PASS: all gamma sweeps completed"
Write-Host "RESULTS=$root"
