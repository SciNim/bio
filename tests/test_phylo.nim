import options
import os
import sugar
import streams
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

  test "Check node echoing, for debug purposes":
    let tree = parse(rootOnLeaf)

    check $tree.nodes[0] == "B:0.2"
    check $tree.nodes[3] == "(C:0.3,D:0.4)E:0.5"

    let tree2 = parse(namedAll)
    check $tree2.nodes[0] == "A"
    check $tree2.nodes[4] == "(C,D)E"
    check $tree2.nodes[5] == "(A,B,(C,D)E)F"

  test "Tree echoing":
    let tree = parse(rootOnLeaf)
    check $tree == "((B:0.2,(C:0.3,D:0.4)E:0.5)F:0.1)A;"

    let tree2 = parse(noNames)
    check $tree2 == "(,,(,));"

    let tree3 = parse("")
    check $tree3 == ""

  test "Find nodes by label":
    let tree = parse(rootOnLeaf)
    check $tree["E"] == "(C:0.3,D:0.4)E:0.5"

    let tree2 = parse(namedAll)
    check $tree2["E"] == "(C,D)E"
    expect KeyError:
      discard $tree2["None"]

suite "Test Newick traverses":
  setup:
    let tree = parse("((A1,A2)B,(C,D)E)F;")

  test "Breath First traverse":
    let nodesA = collect(newSeq):
      for node in traverseBF(tree["A1"]):  # No children
        node.label
    check nodesA == @["A1"]

    let nodesC = collect(newSeq):
      for node in traverseBF(tree["E"]):  # Two children
        node.label
    check nodesC == @["E", "C", "D"]

    let nodesE = collect(newSeq):
      for node in traverseBF(tree["F"]):  # Five children
        node.label
    check nodesE == @["F", "B", "E", "A1", "A2", "C", "D"]

  test "Depth First traverse":
    let nodesA = collect(newSeq):
      for node in traverseDF(tree["A1"]):  # No children
        node.label
    check nodesA == @["A1"]

    let nodesC = collect(newSeq):
      for node in traverseDF(tree["E"]):  # Node E, two children
        node.label
    check nodesC == @["E", "D", "C"]

    let nodesE = collect(newSeq):
      for node in traverseDF(tree["F"]):  # Node E, two children
        node.label
    check nodesE == @["F", "E", "D", "C", "B", "A2", "A1"]
