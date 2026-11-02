param(
    [string]$Python = "python",
    [string]$CellposePython = (Join-Path $env:USERPROFILE "cellpose_env\Scripts\python.exe"),
    [string]$Ilastik = "C:\Program Files\ilastik-1.4.2\ilastik.exe",
    [string]$IlastikPython = "C:\Program Files\ilastik-1.4.2\python.exe"
)
$ErrorActionPreference = "Stop"
$parent = Split-Path -Parent $PSScriptRoot
$environmentPath = Join-Path $parent "granuscape_env"
$temporaryPath = Join-Path $parent "granuscape_temp"
New-Item -ItemType Directory -Force -Path $temporaryPath | Out-Null
$previousTemp, $previousTmp = $env:TEMP, $env:TMP
$env:TEMP = $temporaryPath
$env:TMP = $temporaryPath
try {
    if (-not (Test-Path (Join-Path $environmentPath "Scripts\python.exe"))) {
        & $Python -m venv $environmentPath
        if ($LASTEXITCODE -ne 0) { throw "Could not create the GranuScape environment" }
    }
    $environmentPython = Join-Path $environmentPath "Scripts\python.exe"
    & $environmentPython -m pip install --no-cache-dir -e "$($PSScriptRoot)[test]"
    if ($LASTEXITCODE -ne 0) { throw "Package installation failed" }
    & $environmentPython -m pytest (Join-Path $PSScriptRoot "tests") -q
    if ($LASTEXITCODE -ne 0) { throw "Automated checks failed" }
    New-Item -ItemType Directory -Force -Path (Join-Path $PSScriptRoot "data\test") | Out-Null
    $configurationPath = Join-Path $PSScriptRoot "test.yaml"
    if (-not (Test-Path $configurationPath)) {
        $configuration = [ordered]@{
            input = "./data/test"
            output = "./results/test_01"
            channels = @{dapi="DAPI"; af488="AF488"; rhoda="Rhoda"}
            backends = @{cellpose_python=$CellposePython; cellpose_model="cpsam_v2"; ilastik=$Ilastik; ilastik_python=$IlastikPython}
            gpu = $true
            keep_work = $true
            density_sigma_um = 12
            qc = @{views=@("A","B","C","D"); include_cell_ids=$true; normalization="clahe"; mcl_color="#C86F80"; mcl_line_width_px=3}
            experiment = @{name="integration_test"}
        }
        $json = $configuration | ConvertTo-Json -Depth 6
        [System.IO.File]::WriteAllText($configurationPath, $json, (New-Object System.Text.UTF8Encoding($false)))
    }
    Write-Host "Setup complete. Copy 3-5 test CZI files into data\test, check test.yaml, then run CHECK_TEST.bat."
} finally {
    $env:TEMP = $previousTemp
    $env:TMP = $previousTmp
}
