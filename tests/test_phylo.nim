import options
import os
import sugar
import streams
import strutils
import unittest

import bio / [phylo / newick]


suite "Test Newick parsing":
  setup:
    let basicTree: string = currentSourcePath().parentDir / "test_files" /
      "newick_tree.nwk"
    let labeledTree: string = currentSourcePath().parentDir / "test_files" /
      "node_labels.nwk"

    ## These trees come from the Wikipedia entry "Newick_format"
    ## All the same tree, but with different values in length and label
    let noNames = "(,,(,));"
    let namedLeafs = "(A,B,(C,D));"
    let namedAll = "(A,B,(C,D)E)F;"
    let distancesButRoot = "(:0.1,:0.2,(:0.3,:0.4):0.5);"
    let distancesAll = "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;"
    let namedLeafsAndDistances = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
    let namedAllAndDistances = "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;"
    let rootOnLeaf = "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"

  test "Basic trees":
    for tree in [noNames, namedLeafs, namedAll, distancesButRoot, distancesAll,
                 namedLeafsAndDistances, namedAllAndDistances, rootOnLeaf]:
      let tree = parse(noNames)
      check len(tree) == 6

  test "Labels in basic trees":
    proc labels(tree: Tree): seq[string] =
      collect(newSeq):
        for node in tree.nodes:
          node.label

    for src in [distancesButRoot, distancesAll]:
      let tree = parse(src)
      let lbls = labels(tree)
      check lbls == @["", "", "", "", "", ""]

    for src in [namedLeafs, namedLeafsAndDistances]:
      let tree = parse(src)
      let lbls = labels(tree)
      check lbls == @["A", "B", "C", "D", "", ""]

    for src in [namedAll, namedAllAndDistances]:
      let tree = parse(src)
      let lbls = labels(tree)
      check lbls == @["A", "B", "C", "D", "E", "F"]

    let tree = parse(rootOnLeaf)
    let lbls = labels(tree)
    check lbls == @["B", "C", "D", "E", "F", "A"]

  test "Lengths in basic trees":
    proc lengths(tree: Tree): seq[Option[float]] =
      collect(newSeq):
        for node in tree.nodes:
          node.length

    var tree = parse(namedLeafs)
    var lens = lengths(tree)
    check lens == @[none(float), none(float), none(float),
                    none(float), none(float), none(float)]

    tree = parse(distancesButRoot)
    lens = lengths(tree)
    check lens == @[option(0.1), option(0.2), option(0.3),
                    option(0.4), option(0.5), none(float)]

    tree = parse(distancesAll)
    lens = lengths(tree)
    check lens == @[option(0.1), option(0.2), option(0.3),
                    option(0.4), option(0.5), option(0.0)]

  test "The double single quote in a quoted label":
    let weirdLabel = "(,,(,'A label with ''subquotes'''));"
    let tree = parse(weirdLabel)

    check tree.nodes[3].label == "A label with 'subquotes'"

  test "Parse a File Stream":
    proc labels(tree: Tree): seq[string] =
      collect(newSeq):
        for node in tree.nodes:
          node.label

    proc comments(tree: Tree): seq[string] =
      collect(newSeq):
        for node in tree.nodes:
          node.comment

    var src = newFileStream(basicTree)
    defer: src.close

    let tree = parse(src)
    let lbls = labels(tree)
    let cmts = comments(tree)
    check tree.len == 72

    check "Rangifer tarandus" in lbls  # This label is broken by a \n
    check "Equus ferus caballus" notin lbls  # This is a label-like but in a
                                             # comment
    check "brown bear" in cmts
    check "wild horse; also 'Equus ferus caballus'" in cmts
