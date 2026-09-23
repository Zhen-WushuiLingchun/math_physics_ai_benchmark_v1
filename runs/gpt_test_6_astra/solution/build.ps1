$ErrorActionPreference = 'Stop'
$workspacePath = Split-Path -Parent $PSScriptRoot
$sourcePath = Join-Path $PSScriptRoot 'q1_solution.tex'
$outputPath = Join-Path $workspacePath 'output\pdf'
New-Item -ItemType Directory -Force -Path $outputPath | Out-Null
Push-Location $workspacePath
try {
    python (Join-Path $PSScriptRoot 'checks\verify_q1.py')
    if ($LASTEXITCODE -ne 0) { throw 'Exact verification failed.' }
    for ($pass = 1; $pass -le 3; $pass++) {
        & xelatex '-interaction=nonstopmode' '-halt-on-error' "-output-directory=$outputPath" $sourcePath
        if ($LASTEXITCODE -ne 0) { throw "XeLaTeX pass $pass failed." }
    }
} finally { Pop-Location }
