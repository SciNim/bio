import os
import strformat
# Package

version       = "0.2.9"
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
  let media = ["logo.svg", "PhredVsError.svg"]
  for img in media:
    if fileExists(&"media/{img}"):
      cpFile(&"media/{img}", &"{docDir}/{img}")
    if fileExists(&"src/media/{img}"):
      cpFile(&"src/media/{img}", &"{docDir}/{img}")

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

  for nim_src in ["sequences", "io", "fasta", "fastq", "entrez"]:
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
