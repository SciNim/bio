import os
import strformat
# Package

version       = "0.2.7"
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

const docDir = "public"
proc docSetup =
  mkdir(docDir)
  if fileExists("media/logo.svg"):
    cpFile("media/logo.svg", &"{docDir}/logo.svg")
  if fileExists("src/media/logo.svg"):
    cpFile("src/media/logo.svg", &"{docDir}/logo.svg")

proc preDocs =
  exec("rm *idx -f")
  putEnv("author", author)
  putEnv("nimversion", system.NimVersion)
  putEnv("libversion", version)

proc buildDocs(rst, src: string) =
  for rst_src in ["index.rst", "tutorial.rst", "recipes.rst"]:
    # TODO: Add a "if newer" here.
    let outHtml = rst_src.replace(".rst", ".html")
    selfExec(&"-o:{outHtml} rst2html {rst}/{rst_src}")

  for nim_src in ["sequences", "fasta", "entrez"]:
    selfExec(&"doc --project --index:on --outdir:. {src}/{nim_src}.nim")

proc postDocs() =
  selfExec("buildIndex -o:theindex.html .")
  # Download the JS only if it's newest than the local one.
  exec("curl -LO -z dochack.js " & "https://nim-lang.github.io/Nim/dochack.js")

task docs, &"Deploy doc html + search index to {docDir}/ directory":
  docSetup()

  withDir docDir:
    preDocs()
    buildDocs(rst = "../docs", src = "../bio")
    postDocs()

task repodocs, "Deploy docs, but from repo":
  docSetup()

  withDir docDir:
    preDocs()
    buildDocs(rst = "../src/docs", src = "../src/bio")
    postDocs()
