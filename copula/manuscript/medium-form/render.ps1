$ErrorActionPreference = 'Stop'
$rscript = 'C:\Program Files\R\R-4.5.3\bin\Rscript.exe'
$env:LANG = ''
$env:LANGUAGE = ''
$env:LC_ALL = ''
$env:RSTUDIO_PANDOC = 'C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools'

Push-Location $PSScriptRoot
try {
  & $rscript --vanilla -e "rmarkdown::render('article-medium.Rmd', clean=TRUE, quiet=FALSE)"
  if ($LASTEXITCODE -ne 0) { throw 'Medium-form article render failed' }
} finally {
  Pop-Location
}
Write-Host 'Rendered medium-form/article-medium.pdf'
