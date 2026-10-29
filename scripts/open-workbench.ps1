param([switch]$NoBrowser)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$baseUrl = 'http://127.0.0.1:8767'
$savedJob = 'acceptance-1000-20261027'
$pageUrl = if (Test-Path -LiteralPath (Join-Path $projectRoot "web\jobs\$savedJob\meta.json")) {
    "$baseUrl/?job=$savedJob"
} else {
    "$baseUrl/"
}
$mutex = New-Object System.Threading.Mutex($false, 'Local\ReaxFF_Workbench_8767')
$lockTaken = $false

function Test-Workbench {
    try {
        $health = Invoke-RestMethod -Uri "$baseUrl/healthz" -TimeoutSec 2
        return ($health.status -eq 'ok' -and $health.streaming_analysis -eq $true -and $health.local_path_enabled -eq $true)
    } catch {
        return $false
    }
}

try {
    try { $lockTaken = $mutex.WaitOne(30000) }
    catch [System.Threading.AbandonedMutexException] { $lockTaken = $true }
    if (-not $lockTaken) { throw 'Another launch is still in progress. Please try again shortly.' }

    if (-not (Test-Workbench)) {
        $listener = Get-NetTCPConnection -State Listen -LocalPort 8767 -ErrorAction SilentlyContinue
        if ($listener) { throw 'Port 8767 is occupied but the workbench is not responding. No process was stopped. Please retry or check the service.' }

        $pythonCandidates = @(
            (Join-Path $projectRoot '.venv\Scripts\python.exe'),
            'D:\Users\MECHREVO\miniconda3\python.exe'
        )
        $pythonPath = $pythonCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
        if (-not $pythonPath) {
            $pythonCommand = Get-Command python.exe -ErrorAction SilentlyContinue
            if ($pythonCommand) { $pythonPath = $pythonCommand.Source }
        }
        if (-not $pythonPath) { throw 'Python was not found. Check the Python installation path in scripts\open-workbench.ps1.' }

        $runtimeDir = Join-Path $projectRoot 'web\.runtime'
        New-Item -ItemType Directory -Path $runtimeDir -Force | Out-Null
        $stamp = Get-Date -Format 'yyyyMMdd-HHmmss-ffff'
        $stdoutPath = Join-Path $runtimeDir "workbench-$stamp.stdout.log"
        $stderrPath = Join-Path $runtimeDir "workbench-$stamp.stderr.log"
        $previousPythonPath = $env:PYTHONPATH
        try {
            $env:PYTHONPATH = "$projectRoot;$(Join-Path $projectRoot 'web')"
            if ($previousPythonPath) { $env:PYTHONPATH += ";$previousPythonPath" }
            $server = Start-Process -FilePath $pythonPath -ArgumentList @(
                '-m', 'uvicorn', 'app.main:app', '--host', '127.0.0.1', '--port', '8767'
            ) -WorkingDirectory $projectRoot -WindowStyle Hidden -PassThru `
                -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
        } finally {
            $env:PYTHONPATH = $previousPythonPath
        }

        $ready = $false
        $deadline = [DateTime]::UtcNow.AddSeconds(60)
        while ([DateTime]::UtcNow -lt $deadline) {
            if (Test-Workbench) { $ready = $true; break }
            $server.Refresh()
            if ($server.HasExited) { throw "The service exited during startup. See $stderrPath" }
            Start-Sleep -Milliseconds 500
        }
        if (-not $ready) { throw "The service is still starting. Retry shortly or check $stderrPath" }
        Write-Host "Workbench started (PID $($server.Id))."
    } else {
        Write-Host 'Workbench is already running.'
    }

    if (-not $NoBrowser) { Start-Process -FilePath $pageUrl }
    Write-Host $pageUrl
} catch {
    Write-Host "Could not open the workbench: $($_.Exception.Message)" -ForegroundColor Red
    exit 1
} finally {
    if ($lockTaken) { $mutex.ReleaseMutex() }
    $mutex.Dispose()
}
