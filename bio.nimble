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

  mkdir("htmldocs")

  withDir "htmldocs":
    exec("rm *idx -f")
    putEnv("author", author)
    putEnv("nimversion", system.NimVersion)
    putEnv("libversion", version)
    selfExec(&"rst2html ../docs/index.rst")
    selfExec("doc --project --index:on ../src/bio/sequences.nim")
    selfExec("doc --project --index:on ../src/bio/fasta.nim")
    selfExec("buildIndex -o:theindex.html .")
    # Download the JS only if it's newest than the local one.
    exec("curl -LO -z dochack.js " & docHackJsSource)
