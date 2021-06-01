import os
import strformat
# Package

version       = "0.2.6"
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

proc docSetup =
  mkdir("htmldocs")
  cpFile("media/logo.svg", "htmldocs/logo.svg")

proc preDocs =
  exec("rm *idx -f")
  putEnv("author", author)
  putEnv("nimversion", system.NimVersion)
  putEnv("libversion", version)

proc buildDocs(rst, src: string) =
  for rst_src in ["index.rst", "tutorial.rst", "recipes.rst"]:
    # TODO: Add a "if newer" here.
    selfExec(&"rst2html {rst}/{rst_src}")

  for nim_src in ["sequences", "fasta", "entrez"]:
    selfExec(&"doc --project --index:on {src}/{nim_src}.nim")

  for outputFile in listFiles(rst & "/htmldocs"):
    mvfile(outputFile, &"{splitPath(outputFile)[1]}")

  for outputFile in listFiles(src & "/htmldocs"):
    mvfile(outputFile, &"{splitPath(outputFile)[1]}")

  rmDir(&"{rst}/htmldocs")
  rmDir(&"{src}/htmldocs")


proc postDocs() =
  selfExec("buildIndex -o:theindex.html .")
  # Download the JS only if it's newest than the local one.
  exec("curl -LO -z dochack.js " & "https://nim-lang.github.io/Nim/dochack.js")

task docs, "Deploy doc html + search index to public/ directory":
  docSetup()

  withDir "htmldocs":
    preDocs()
    buildDocs(rst = "../docs", src = "../bio")
    postDocs()

task repodocs, "Deploy docs, but from repo":
  docSetup()

  withDir "htmldocs":
    preDocs()
    buildDocs(rst = "../src/docs", src = "../src/bio")
    cpFile("../media/logo.svg", "logo.svg")
    postDocs()
