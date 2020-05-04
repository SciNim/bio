import strformat
# Package

version       = "0.1.0"
author        = "Xabier Bello"
description   = "A set of tools to do biological operations."
license       = "MIT"
srcDir        = "src"

# Dependencies

requires "nim >= 1.2.0"

# Tasks
task test, "Full test suite":
  exec "testament p tests"
  #exec "testament html"
  #exec "firefox testresults.html"

task docs, "Deploy doc html + search index to public/ directory":
  let
    docHackJsSource = "https://nim-lang.github.io/Nim/dochack.js"
    envs = &"--putenv:nimversion={system.NimVersion} " &
           &"--putenv:libversion={version}"

  withDir "htmldocs":
    exec("rm *idx -f")
    exec(&"nim rst2html {envs} ../docs/index.rst")
    exec("nim doc --project --index:on ../src/bio.nim")
    exec("nim doc --index:on ../src/bio/fasta.nim")
    exec("nim buildIndex -o:theindex.html .")
    # Download the JS only if it's newest than the local one.
    exec("curl -LO -z dochack.js " & docHackJsSource)
