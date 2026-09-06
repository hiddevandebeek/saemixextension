param(
  [Parameter(Mandatory=$true)][string]$Tag,
  [ValidateRange(0,60000)][int]$DelayMilliseconds = 0,
  [ValidateRange(0,255)][int]$ExitCode = 0,
  [Parameter(Mandatory=$true)][string]$ResultPath,
  [string]$ResumeToken = ""
)

$ErrorActionPreference = "Stop"

Write-Output ("STDOUT|{0}" -f $Tag)
[Console]::Error.WriteLine(("STDERR|{0}" -f $Tag))

if ($DelayMilliseconds -gt 0) {
  Start-Sleep -Milliseconds $DelayMilliseconds
}

$resultFull = [System.IO.Path]::GetFullPath($ResultPath)
$resultDirectory = [System.IO.Path]::GetDirectoryName($resultFull)
if ([string]::IsNullOrWhiteSpace($resultDirectory)) {
  throw "ResultPath must have a parent directory"
}
[System.IO.Directory]::CreateDirectory($resultDirectory) | Out-Null
$utf8NoBom = New-Object System.Text.UTF8Encoding($false)

if ([System.IO.File]::Exists($resultFull)) {
  if ([string]::IsNullOrEmpty($ResumeToken)) {
    [Console]::Error.WriteLine(("INCOMPATIBLE|{0}|missing-resume-token" -f $Tag))
    exit 86
  }
  $existing = [System.IO.File]::ReadAllText($resultFull, $utf8NoBom)
  if ($existing -ceq $ResumeToken) {
    Write-Output ("SKIP|{0}" -f $Tag)
    exit 0
  }
  [Console]::Error.WriteLine(("INCOMPATIBLE|{0}|token-mismatch" -f $Tag))
  exit 86
}

$token = if ([string]::IsNullOrEmpty($ResumeToken)) {
  "RESULT|{0}" -f $Tag
} else {
  $ResumeToken
}
$temporary = "{0}.tmp.{1}.{2}" -f $resultFull,$PID,[guid]::NewGuid().ToString("N")
try {
  [System.IO.File]::WriteAllText($temporary, $token, $utf8NoBom)
  [System.IO.File]::Move($temporary, $resultFull)
} finally {
  if ([System.IO.File]::Exists($temporary)) {
    [System.IO.File]::Delete($temporary)
  }
}

Write-Output ("WROTE|{0}" -f $Tag)
exit $ExitCode
