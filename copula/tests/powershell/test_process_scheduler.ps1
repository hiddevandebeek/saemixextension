param(
  [string]$ShellPath = "",
  [switch]$BehavioralOnly,
  [switch]$StaticOnly
)

$ErrorActionPreference = "Stop"
if ($BehavioralOnly -and $StaticOnly) {
  throw "BehavioralOnly and StaticOnly are mutually exclusive"
}

$runBehavioral = -not $StaticOnly
$runStatic = -not $BehavioralOnly
$script:Failures = New-Object 'System.Collections.Generic.List[object]'
$fakeChild = (Resolve-Path (Join-Path $PSScriptRoot "fake_child.ps1")).Path
$repo = (Resolve-Path (Join-Path $PSScriptRoot "..\..\..")).Path

if ([string]::IsNullOrWhiteSpace($ShellPath)) {
  $ShellPath = (Get-Process -Id $PID).Path
}
$ShellPath = [System.IO.Path]::GetFullPath($ShellPath)
if (-not [System.IO.File]::Exists($ShellPath)) {
  throw "ShellPath does not exist: $ShellPath"
}

function Add-Failure {
  param([string]$Case,[string]$Assertion,[string]$Detail)
  [void]$script:Failures.Add([pscustomobject]@{
    Case=$Case; Assertion=$Assertion; Detail=$Detail
  })
  Write-Output ("FAIL|{0}|{1}|{2}" -f $Case,$Assertion,$Detail)
}

function Assert-Condition {
  param([bool]$Condition,[string]$Case,[string]$Assertion,[string]$Detail="")
  if ($Condition) {
    Write-Output ("PASS|{0}|{1}" -f $Case,$Assertion)
  } else {
    Add-Failure -Case $Case -Assertion $Assertion -Detail $Detail
  }
}

function Quote-PowerShellLiteral {
  param([AllowEmptyString()][string]$Value)
  return "'" + $Value.Replace("'", "''") + "'"
}

function New-PlanItem {
  param(
    [string]$Tag,
    [int]$DelayMilliseconds,
    [int]$ExitCode,
    [string]$ResultPath,
    [string]$ResumeToken = ""
  )
  return [pscustomobject]@{
    Tag=$Tag
    DelayMilliseconds=$DelayMilliseconds
    ExitCode=$ExitCode
    ResultPath=$ResultPath
    ResumeToken=$ResumeToken
  }
}

function New-ChildArguments {
  param([object]$Item)
  $command = "& {0} -Tag {1} -DelayMilliseconds {2} -ExitCode {3} -ResultPath {4} -ResumeToken {5}" -f
    (Quote-PowerShellLiteral $fakeChild),
    (Quote-PowerShellLiteral ([string]$Item.Tag)),
    ([int]$Item.DelayMilliseconds),
    ([int]$Item.ExitCode),
    (Quote-PowerShellLiteral ([string]$Item.ResultPath)),
    (Quote-PowerShellLiteral ([string]$Item.ResumeToken))
  $encoded = [Convert]::ToBase64String([Text.Encoding]::Unicode.GetBytes($command))
  return @("-NoLogo","-NoProfile","-NonInteractive","-ExecutionPolicy","Bypass",
    "-EncodedCommand",$encoded)
}

function Invoke-SchedulerFixture {
  param(
    [object[]]$Plan,
    [int]$Workers,
    [string]$CaseRoot,
    [int]$PollMilliseconds = 15
  )
  if ($Workers -lt 1) { throw "Workers must be positive" }
  [System.IO.Directory]::CreateDirectory($CaseRoot) | Out-Null
  $logRoot = Join-Path $CaseRoot "logs"
  [System.IO.Directory]::CreateDirectory($logRoot) | Out-Null

  $queue = New-Object 'System.Collections.Generic.Queue[object]'
  foreach ($item in @($Plan)) { $queue.Enqueue($item) }
  $active = @()
  $events = New-Object 'System.Collections.Generic.List[object]'
  $statuses = New-Object 'System.Collections.Generic.List[object]'
  $processObjects = New-Object 'System.Collections.Generic.List[object]'
  $allPids = New-Object 'System.Collections.Generic.List[int]'
  $clock = [System.Diagnostics.Stopwatch]::StartNew()
  $sequence = 0
  $pollCycle = 0
  $maxActive = 0
  $fixtureResult = $null

  try {
    while ($queue.Count -gt 0 -or $active.Count -gt 0) {
      while ($queue.Count -gt 0 -and $active.Count -lt $Workers) {
        $item = $queue.Dequeue()
        $stdout = Join-Path $logRoot (([string]$item.Tag) + ".out.log")
        $stderr = Join-Path $logRoot (([string]$item.Tag) + ".err.log")
        $arguments = New-ChildArguments -Item $item
        $beforeTags = @($active | ForEach-Object { $_.Tag })
        $process = Start-Process -FilePath $ShellPath -ArgumentList $arguments `
          -PassThru -WindowStyle Hidden -WorkingDirectory $CaseRoot `
          -RedirectStandardOutput $stdout -RedirectStandardError $stderr
        try {
          $null = $process.Handle
        } catch {
          try {
            if (-not $process.HasExited) { $process.Kill() }
            $process.WaitForExit()
          } finally {
            $process.Dispose()
          }
          throw "Could not retain process handle for $($item.Tag), pid=$($process.Id): $($_.Exception.Message)"
        }

        [void]$processObjects.Add($process)
        [void]$allPids.Add([int]$process.Id)
        $active += [pscustomobject]@{
          Process=$process
          Item=$item
          Tag=[string]$item.Tag
          StdoutPath=$stdout
          StderrPath=$stderr
        }
        if ($active.Count -gt $maxActive) { $maxActive = $active.Count }
        $sequence++
        [void]$events.Add([pscustomobject]@{
          Sequence=$sequence
          Kind="LAUNCH"
          Tag=[string]$item.Tag
          Pid=[int]$process.Id
          ElapsedMilliseconds=[int64]$clock.ElapsedMilliseconds
          PollCycle=$pollCycle
          ActiveCount=$active.Count
          ActiveTagsBefore=$beforeTags
        })
      }

      if ($active.Count -gt 0) {
        Start-Sleep -Milliseconds $PollMilliseconds
      }
      $pollCycle++
      $still = @()
      foreach ($job in $active) {
        if ($job.Process.HasExited) {
          $job.Process.WaitForExit()
          $exitValue = $job.Process.ExitCode
          if ($null -eq $exitValue) {
            throw "Null exit code for $($job.Tag), pid=$($job.Process.Id)"
          }
          [int]$exitCode = $exitValue
          $status = [pscustomobject]@{
            Tag=$job.Tag
            Pid=[int]$job.Process.Id
            ExitCode=$exitCode
            Classification=if ($exitCode -eq 0) { "DONE" } else { "FAIL" }
            StdoutPath=$job.StdoutPath
            StderrPath=$job.StderrPath
            ResultPath=[string]$job.Item.ResultPath
            Disposed=$false
          }
          try {
            [void]$statuses.Add($status)
            $sequence++
            [void]$events.Add([pscustomobject]@{
              Sequence=$sequence
              Kind="COMPLETE"
              Tag=$job.Tag
              Pid=[int]$job.Process.Id
              ExitCode=$exitCode
              ElapsedMilliseconds=[int64]$clock.ElapsedMilliseconds
              PollCycle=$pollCycle
              ActiveCount=$active.Count
              ActiveTagsBefore=@()
            })
          } finally {
            $job.Process.Dispose()
            $status.Disposed = $true
          }
        } else {
          $still += $job
        }
      }
      $active = $still
    }

    $aggregate = if (@($statuses | Where-Object { $_.ExitCode -ne 0 }).Count -gt 0) { 1 } else { 0 }
    $fixtureResult = [pscustomobject]@{
      Statuses=@($statuses)
      Events=@($events)
      MaximumActive=$maxActive
      AggregateExitCode=$aggregate
      ProcessObjects=@($processObjects)
      Pids=@($allPids)
      PollCycles=$pollCycle
    }
  } finally {
    foreach ($job in @($active)) {
      try {
        if (-not $job.Process.HasExited) { $job.Process.Kill() }
        $job.Process.WaitForExit()
      } catch {
      } finally {
        $job.Process.Dispose()
      }
    }
    $clock.Stop()
  }
  return $fixtureResult
}

function Assert-LogIsolation {
  param([string]$Case,[object]$Result,[string[]]$Tags)
  foreach ($status in @($Result.Statuses)) {
    $stdoutLines = @(Get-Content -LiteralPath $status.StdoutPath -ErrorAction Stop)
    $stderrLines = @(Get-Content -LiteralPath $status.StderrPath -ErrorAction Stop)
    Assert-Condition -Condition ($stdoutLines -contains ("STDOUT|{0}" -f $status.Tag)) `
      -Case $Case -Assertion ("stdout-complete-{0}" -f $status.Tag)
    Assert-Condition -Condition ($stderrLines -contains ("STDERR|{0}" -f $status.Tag)) `
      -Case $Case -Assertion ("stderr-complete-{0}" -f $status.Tag)
    foreach ($other in $Tags) {
      if ($other -ceq $status.Tag) { continue }
      $foreignOut = @($stdoutLines | Where-Object { $_ -match ("\|{0}(\||$)" -f [regex]::Escape($other)) })
      $foreignErr = @($stderrLines | Where-Object { $_ -match ("\|{0}(\||$)" -f [regex]::Escape($other)) })
      Assert-Condition -Condition (($foreignOut.Count + $foreignErr.Count) -eq 0) `
        -Case $Case -Assertion ("log-isolation-{0}-not-{1}" -f $status.Tag,$other)
    }
  }
}

function Test-DisposedProcess {
  param([System.Diagnostics.Process]$Process)
  try {
    $null = $Process.Handle
    return $false
  } catch [System.ObjectDisposedException] {
    return $true
  } catch [System.InvalidOperationException] {
    return $true
  }
}

function Invoke-BehavioralTests {
  param([string]$TemporaryRoot)

  $case = "fast-success"
  $root = Join-Path $TemporaryRoot $case
  $plan = @(foreach ($i in 1..25) {
    $tag = "fast_{0:D2}" -f $i
    New-PlanItem -Tag $tag -DelayMilliseconds 0 -ExitCode 0 `
      -ResultPath (Join-Path $root ("results\{0}.token" -f $tag))
  })
  try {
    $result = Invoke-SchedulerFixture -Plan $plan -Workers 6 -CaseRoot $root
    Assert-Condition -Condition (@($result.Statuses).Count -eq 25) -Case $case -Assertion "all-25-drained"
    Assert-Condition -Condition (@($result.Statuses | Where-Object { $null -eq $_.ExitCode }).Count -eq 0) `
      -Case $case -Assertion "no-null-exit-code"
    Assert-Condition -Condition (@($result.Statuses | Where-Object { $_.ExitCode -ne 0 }).Count -eq 0) `
      -Case $case -Assertion "all-exit-zero"
    Assert-Condition -Condition (@($result.Statuses | Where-Object { $_.Classification -ne "DONE" }).Count -eq 0) `
      -Case $case -Assertion "all-classified-done"
    Assert-Condition -Condition ($result.AggregateExitCode -eq 0) -Case $case -Assertion "aggregate-zero"
    Assert-LogIsolation -Case $case -Result $result -Tags @($plan.Tag)
  } catch {
    Add-Failure -Case $case -Assertion "case-execution" -Detail $_.Exception.Message
  }

  $case = "real-failure"
  $root = Join-Path $TemporaryRoot $case
  $codes = @(0,7,0,23,0)
  $plan = @(for ($i=0; $i -lt $codes.Count; $i++) {
    $tag = "failure_{0:D2}" -f ($i + 1)
    New-PlanItem -Tag $tag -DelayMilliseconds 0 -ExitCode $codes[$i] `
      -ResultPath (Join-Path $root ("results\{0}.token" -f $tag))
  })
  try {
    $result = Invoke-SchedulerFixture -Plan $plan -Workers 3 -CaseRoot $root
    $failed = @($result.Statuses | Where-Object { $_.Classification -eq "FAIL" })
    Assert-Condition -Condition ($failed.Count -eq 2) -Case $case -Assertion "exactly-two-failures"
    Assert-Condition -Condition ((@($failed.ExitCode | Sort-Object) -join ",") -ceq "7,23") `
      -Case $case -Assertion "exact-nonzero-codes" -Detail (@($failed.ExitCode) -join ",")
    Assert-Condition -Condition (@($result.Statuses).Count -eq $plan.Count) -Case $case -Assertion "all-children-drained"
    Assert-Condition -Condition ($result.AggregateExitCode -eq 1) -Case $case -Assertion "aggregate-one"
    Assert-LogIsolation -Case $case -Result $result -Tags @($plan.Tag)
  } catch {
    Add-Failure -Case $case -Assertion "case-execution" -Detail $_.Exception.Message
  }

  $case = "scheduling-refill"
  $root = Join-Path $TemporaryRoot $case
  $plan = @(
    (New-PlanItem -Tag "slow" -DelayMilliseconds 800 -ExitCode 0 -ResultPath (Join-Path $root "results\slow.token")),
    (New-PlanItem -Tag "fast" -DelayMilliseconds 0 -ExitCode 0 -ResultPath (Join-Path $root "results\fast.token")),
    (New-PlanItem -Tag "fast_fail" -DelayMilliseconds 0 -ExitCode 7 -ResultPath (Join-Path $root "results\fast_fail.token")),
    (New-PlanItem -Tag "short_1" -DelayMilliseconds 40 -ExitCode 0 -ResultPath (Join-Path $root "results\short_1.token")),
    (New-PlanItem -Tag "short_2" -DelayMilliseconds 50 -ExitCode 0 -ResultPath (Join-Path $root "results\short_2.token"))
  )
  try {
    $result = Invoke-SchedulerFixture -Plan $plan -Workers 2 -CaseRoot $root
    $launches = @($result.Events | Where-Object { $_.Kind -eq "LAUNCH" } | Sort-Object Sequence)
    $completions = @($result.Events | Where-Object { $_.Kind -eq "COMPLETE" })
    Assert-Condition -Condition ((@($launches.Tag) -join "|") -ceq (@($plan.Tag) -join "|")) `
      -Case $case -Assertion "fifo-launch-order"
    Assert-Condition -Condition ($result.MaximumActive -eq 2) -Case $case -Assertion "maximum-active-exactly-two"
    Assert-Condition -Condition (@($launches | Where-Object { $_.ActiveCount -gt 2 }).Count -eq 0) `
      -Case $case -Assertion "never-exceeds-two"
    $third = $launches[2]
    $slowDone = @($completions | Where-Object { $_.Tag -eq "slow" })[0]
    Assert-Condition -Condition (@($third.ActiveTagsBefore) -contains "slow") `
      -Case $case -Assertion "third-launched-while-slow-tracked"
    Assert-Condition -Condition ($third.ElapsedMilliseconds -lt $slowDone.ElapsedMilliseconds) `
      -Case $case -Assertion "third-launched-before-slow-completed"
    Assert-Condition -Condition ($result.AggregateExitCode -eq 1) -Case $case -Assertion "failure-does-not-stop-drain"
    Assert-LogIsolation -Case $case -Result $result -Tags @($plan.Tag)
  } catch {
    Add-Failure -Case $case -Assertion "case-execution" -Detail $_.Exception.Message
  }

  $case = "resume"
  $root = Join-Path $TemporaryRoot $case
  [System.IO.Directory]::CreateDirectory((Join-Path $root "results")) | Out-Null
  $utf8NoBom = New-Object System.Text.UTF8Encoding($false)
  $compatiblePath = Join-Path $root "results\compatible.token"
  $incompatiblePath = Join-Path $root "results\incompatible.token"
  [System.IO.File]::WriteAllText($compatiblePath, "token-v1", $utf8NoBom)
  [System.IO.File]::WriteAllText($incompatiblePath, "wrong-token", $utf8NoBom)
  $compatibleBefore = [Convert]::ToBase64String([System.IO.File]::ReadAllBytes($compatiblePath))
  $incompatibleBefore = [Convert]::ToBase64String([System.IO.File]::ReadAllBytes($incompatiblePath))
  $plan = @(
    (New-PlanItem -Tag "resume_ok" -DelayMilliseconds 0 -ExitCode 0 `
      -ResultPath $compatiblePath -ResumeToken "token-v1"),
    (New-PlanItem -Tag "resume_bad" -DelayMilliseconds 0 -ExitCode 0 `
      -ResultPath $incompatiblePath -ResumeToken "token-v1"),
    (New-PlanItem -Tag "resume_refill" -DelayMilliseconds 0 -ExitCode 0 `
      -ResultPath (Join-Path $root "results\refill.token") -ResumeToken "token-v1")
  )
  try {
    $result = Invoke-SchedulerFixture -Plan $plan -Workers 1 -CaseRoot $root
    $ok = @($result.Statuses | Where-Object { $_.Tag -eq "resume_ok" })[0]
    $bad = @($result.Statuses | Where-Object { $_.Tag -eq "resume_bad" })[0]
    $okOut = @(Get-Content -LiteralPath $ok.StdoutPath)
    Assert-Condition -Condition ($ok.ExitCode -eq 0 -and $ok.Classification -eq "DONE") `
      -Case $case -Assertion "compatible-token-exits-zero"
    Assert-Condition -Condition ($okOut -contains "SKIP|resume_ok") -Case $case -Assertion "compatible-token-prints-skip"
    Assert-Condition -Condition ($bad.ExitCode -eq 86 -and $bad.Classification -eq "FAIL") `
      -Case $case -Assertion "incompatible-token-nonzero"
    Assert-Condition -Condition (([Convert]::ToBase64String([System.IO.File]::ReadAllBytes($compatiblePath))) -ceq $compatibleBefore) `
      -Case $case -Assertion "compatible-token-byte-identical"
    Assert-Condition -Condition (([Convert]::ToBase64String([System.IO.File]::ReadAllBytes($incompatiblePath))) -ceq $incompatibleBefore) `
      -Case $case -Assertion "incompatible-token-byte-identical"
    $okComplete = @($result.Events | Where-Object { $_.Kind -eq "COMPLETE" -and $_.Tag -eq "resume_ok" })[0]
    $badLaunch = @($result.Events | Where-Object { $_.Kind -eq "LAUNCH" -and $_.Tag -eq "resume_bad" })[0]
    Assert-Condition -Condition ($badLaunch.PollCycle -eq $okComplete.PollCycle) `
      -Case $case -Assertion "refill-at-next-dispatch-opportunity"
    Assert-Condition -Condition (@($result.Statuses).Count -eq 3) -Case $case -Assertion "resume-queue-drained"
    Assert-LogIsolation -Case $case -Result $result -Tags @($plan.Tag)
  } catch {
    Add-Failure -Case $case -Assertion "case-execution" -Detail $_.Exception.Message
  }

  $case = "repeated-cleanup"
  try {
    foreach ($round in 1..5) {
      $root = Join-Path $TemporaryRoot ("cleanup_{0:D2}" -f $round)
      $plan = @(foreach ($i in 1..3) {
        $tag = "cleanup_{0:D2}_{1:D2}" -f $round,$i
        New-PlanItem -Tag $tag -DelayMilliseconds (($i - 1) * 10) -ExitCode 0 `
          -ResultPath (Join-Path $root ("results\{0}.token" -f $tag))
      })
      $result = Invoke-SchedulerFixture -Plan $plan -Workers 2 -CaseRoot $root
      Assert-Condition -Condition (@($result.Statuses).Count -eq 3) -Case $case `
        -Assertion ("round-{0}-drained" -f $round)
      Assert-Condition -Condition (@($result.Statuses | Where-Object { -not $_.Disposed }).Count -eq 0) `
        -Case $case -Assertion ("round-{0}-status-disposed" -f $round)
      $undisposed = @($result.ProcessObjects | Where-Object { -not (Test-DisposedProcess -Process $_) })
      Assert-Condition -Condition ($undisposed.Count -eq 0) -Case $case `
        -Assertion ("round-{0}-objects-disposed" -f $round)
      $live = @($result.Pids | Where-Object { $null -ne (Get-Process -Id $_ -ErrorAction SilentlyContinue) })
      Assert-Condition -Condition ($live.Count -eq 0) -Case $case `
        -Assertion ("round-{0}-no-live-pids" -f $round) -Detail ($live -join ",")
    }
  } catch {
    Add-Failure -Case $case -Assertion "case-execution" -Detail $_.Exception.Message
  }
}

function Invoke-ProductionShapeTests {
  $launchers = @(
    "copula\study_publication\orchestration\run_primary_balanced.ps1",
    "copula\study_publication\run_manifest.ps1",
    "copula\study_publication\run_calibration.ps1",
    "copula\study_simple\run_manifest.ps1"
  )
  foreach ($relative in $launchers) {
    $case = "production-shape:$relative"
    $path = Join-Path $repo $relative
    $text = [System.IO.File]::ReadAllText($path)
    $starts = [regex]::Matches($text, '(?im)^\s*\$process\s*=\s*Start-Process\b')
    Assert-Condition -Condition ($starts.Count -gt 0) -Case $case -Assertion "has-start-process"
    for ($i=0; $i -lt $starts.Count; $i++) {
      $start = $starts[$i]
      $insert = [regex]::Match($text, '(?im)^\s*\$(?:active|running)\s*\+=', $start.Index + $start.Length)
      $hasInsert = $insert.Success
      Assert-Condition -Condition $hasInsert -Case $case -Assertion ("start-{0}-has-tracked-insertion" -f ($i+1))
      if ($hasInsert) {
        $between = $text.Substring($start.Index, $insert.Index + $insert.Length - $start.Index)
        $handle = [regex]::Match($between, '(?im)\$null\s*=\s*\$process\.Handle')
        Assert-Condition -Condition $handle.Success -Case $case `
          -Assertion ("start-{0}-retains-handle-before-insertion" -f ($i+1))
        $guarded = [regex]::IsMatch($between, '(?is)try\s*\{.*?\$null\s*=\s*\$process\.Handle.*?\}\s*catch')
        Assert-Condition -Condition $guarded -Case $case `
          -Assertion ("start-{0}-handle-acquisition-has-failure-path" -f ($i+1))
      }
    }

    $exits = [regex]::Matches($text, '(?i)\$job\.(?:Process|process)\.HasExited')
    Assert-Condition -Condition ($exits.Count -gt 0) -Case $case -Assertion "has-completion-loop"
    for ($i=0; $i -lt $exits.Count; $i++) {
      $exit = $exits[$i]
      $boundary = [regex]::Match($text, '(?im)^\s*\$(?:active|running)\s*=\s*\$still(?:Running)?\b',
        $exit.Index + $exit.Length)
      if (-not $boundary.Success) {
        Add-Failure -Case $case -Assertion ("completion-{0}-has-boundary" -f ($i+1)) `
          -Detail "could not delimit completion block"
        continue
      }
      $block = $text.Substring($exit.Index, $boundary.Index - $exit.Index)
      $wait = [regex]::Match($block, '(?i)\$job\.(?:Process|process)\.WaitForExit\s*\(')
      $capture = [regex]::Match($block,
        '(?i)\[int\]\s*\$exitCode\s*=\s*\$job\.(?:Process|process)\.ExitCode')
      $dispose = [regex]::Match($block, '(?i)\$job\.(?:Process|process)\.Dispose\s*\(')
      Assert-Condition -Condition $wait.Success -Case $case `
        -Assertion ("completion-{0}-waits" -f ($i+1))
      Assert-Condition -Condition $capture.Success -Case $case `
        -Assertion ("completion-{0}-typed-exit-capture" -f ($i+1))
      Assert-Condition -Condition $dispose.Success -Case $case `
        -Assertion ("completion-{0}-disposes" -f ($i+1))
      $ordered = $wait.Success -and $capture.Success -and $dispose.Success -and
        $wait.Index -lt $capture.Index -and $capture.Index -lt $dispose.Index
      Assert-Condition -Condition $ordered -Case $case `
        -Assertion ("completion-{0}-wait-capture-dispose-order" -f ($i+1))
      $exitReads = [regex]::Matches($block, '(?i)\$job\.(?:Process|process)\.ExitCode').Count
      Assert-Condition -Condition ($exitReads -eq 1) -Case $case `
        -Assertion ("completion-{0}-single-exit-read" -f ($i+1)) -Detail ("reads={0}" -f $exitReads)
    }
    Assert-Condition -Condition (-not [regex]::IsMatch($text, '(?i)\.Refresh\s*\(')) `
      -Case $case -Assertion "no-refresh-workaround"
  }
}

$temporaryRoot = $null
try {
  if ($runBehavioral) {
    $temporaryRoot = Join-Path ([System.IO.Path]::GetTempPath()) `
      ("saemix-copula-scheduler-{0}" -f [guid]::NewGuid().ToString("N"))
    [System.IO.Directory]::CreateDirectory($temporaryRoot) | Out-Null
    Invoke-BehavioralTests -TemporaryRoot $temporaryRoot
  }
  if ($runStatic) {
    Invoke-ProductionShapeTests
  }
} finally {
  if (-not [string]::IsNullOrWhiteSpace($temporaryRoot) -and
      [System.IO.Directory]::Exists($temporaryRoot)) {
    $resolvedTarget = [System.IO.Path]::GetFullPath($temporaryRoot)
    $resolvedBase = [System.IO.Path]::GetFullPath([System.IO.Path]::GetTempPath()).TrimEnd('\','/') +
      [System.IO.Path]::DirectorySeparatorChar
    if (-not $resolvedTarget.StartsWith($resolvedBase, [System.StringComparison]::OrdinalIgnoreCase)) {
      throw "Refusing to remove temporary path outside the system temp directory: $resolvedTarget"
    }
    Remove-Item -LiteralPath $resolvedTarget -Recurse -Force
  }
}

Write-Output ("SUMMARY|shell={0}|behavioral={1}|static={2}|failures={3}" -f
  $ShellPath,$runBehavioral,$runStatic,$script:Failures.Count)
if ($script:Failures.Count -gt 0) {
  exit 1
}
exit 0
