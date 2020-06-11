import strformat
# Package

version       = "0.2.1"
author        = "Xabier Bello"
description   = "A set of tools to do biological operations."
license       = "MIT"
srcDir        = "src"

# Dependencies

requires "nim >= 1.2.0"

# Tasks
#
task test, "Full test suite":
  exec "testament p tests"
  #exec "testament html"
  #exec "firefox testresults.html"

proc preDocs =
  exec("rm *idx -f")
  putEnv("author", author)
  putEnv("nimversion", system.NimVersion)
  putEnv("libversion", version)

proc buildDocs(rst, src: string) =
  for rst_src in ["index.rst", "tutorial.rst"]:
    selfExec(&"rst2html {rst}/{rst_src}")

  for nim_src in ["sequences.nim", "fasta.nim"]:
    selfExec(&"doc --project --index:on {src}/{nim_src}")

proc postDocs() =
  selfExec("buildIndex -o:theindex.html .")
  # Download the JS only if it's newest than the local one.
  exec("curl -LO -z dochack.js " & "https://nim-lang.github.io/Nim/dochack.js")

task docs, "Deploy doc html + search index to public/ directory":
  let
    docHackJsSource = "https://nim-lang.github.io/Nim/dochack.js"

  mkdir("htmldocs")

  withDir "htmldocs":
    preDocs()
    buildDocs(rst = "../docs", src = "../bio")
    postDocs()

task repodocs, "Deploy docs, but from repo":
  mkdir("htmldocs")

  withDir "htmldocs":
    preDocs()
    buildDocs(rst = "../src/docs", src = "../src/bio")
    postDocs()
