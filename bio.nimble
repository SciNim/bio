# Package

version       = "0.1.0"
author        = "Xabier Bello"
description   = "A new awesome nimble package"
license       = "MIT"
srcDir        = "src"

# Dependencies

requires "nim >= 1.2.0"

# Tasks
task test, "Full test suite":
  exec "testament p tests"
  #exec "testament html"
  #exec "firefox testresults.html"
  #
task buildDocs, "Deploy doc html + search index to public/ directory":
  let
    docHackJsSource = "https://nim-lang.github.io/Nim/dochack.js"
  withDir "htmldocs":
    exec("rm *idx -f")
    exec("nim doc --project --index:on ../src/bio.nim")
    #exec("nim doc --project --index:on ../src/bio/fasta.nim")
    #exec("mv bio.html index.html")
    exec("nim buildIndex -o:theindex.html .")
    exec("curl -LO " & docHackJsSource)
    exec("rm *fasta -f")
