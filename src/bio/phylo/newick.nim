## :Author: |author|
## :Version: |libversion|
##
## Parser for Newick files, as defined here_.
##
## This is how the definitions are interpreted:
##
##   - Unquoted labels may not contain some symbols: they are removed.
##   - Single quotes in quoted labels are represented by two single quotes: TODO
##   - Blanks or tabs may appear anywhere except...: they are removed.
##   - Newlines may appear anywhere except...: they are removed.
##
## .. _here: http://evolution.genetics.washington.edu/phylip/newick_doc.html
##
## .. code-block::
##
##   import bio / [phylo / newick]
##
##   for line in "mytree.nh".lines:
##     let tree: Tree = parse(line)
##     echo tree.len
##
##     for node in tree.nodes:
##       echo node.label, ">", node.length
##

import std / [deques, options, strformat, streams, strutils, sugar]


# TODO: types should go to some phylo base file
type
  NodeObj = object of RootObj
    label*, comment*: string
    length*: Option[float]
    parent*: Option[Node]
    children*: seq[Node]
  Node* = ref object of NodeObj
    ## .. code-block::
    ##
    ##    label*, comment*: string
    ##    length*: Option[float]
    ##    parent*: Option[Node]
    ##    children*: seq[Node]

  TreeObj = object of RootObj
    nodes*: seq[Node]
  Tree* = ref object of TreeObj
    ## .. code-block::
    ##
    ##    nodes*: seq[Node]

proc `$`*(n: Node): string =
  ## Returns the string for the Node, including its children.
  ##
  let ownLabel = &"""{(block:
    var r: string = n.label
    if n.length.isSome: r.add ":" & $get(n.length)
    r)}"""

  if len(n.children) == 0:
    return ownLabel
  else:
    return "(" & join(n.children, ",") & ")" & ownLabel

proc `$`*(t: Tree): string =
  result = $(t.nodes[^1])
  if not result.isEmptyOrWhitespace:
    result.add ";"

proc `[]`*(t: Tree, l: string): Node =
  ## Find a Node in the Tree using the label.
  ##
  runnableExamples:
    let tree: Tree = parse("((A,(B1,B2)B)C,D)")
    doAssert $tree["B"] == "(B1,B2)B"

  for node in t.nodes:
    if node.label == l:
      return node
  raise newException(KeyError, "key not found: " & l)

proc len*(t: Tree): Natural =
  len(t.nodes)

iterator traverseBF*(n: Node): Node =
  ## Traverses (searchs) a tree starting at node `n` and yielding nodes using
  ## `Bread First Search`_.
  ##
  ## .. _`Bread First Search`: https://en.wikipedia.org/wiki/Breadth-first_search
  ##
  ## E.g. Given the tree:
  ##
  ## .. Don't delete the empty line below (it's a \u00a0 char), or the tree
  ##    doesn't render
  ## .. code-block::
  ##
  ##   ??
  ##            ??????A1
  ##          ???A???
  ##          ??? ??????A2
  ##        ???C??? ??????B1
  ##      ????????? ???B???
  ##        ???   ??????B2
  ##        ??????????????????D
  ##
  ## Searching Breath First from node C, it yields A, B, A1, A2, B1 and B2, in
  ## that order.
  runnableExamples:
    let tree: Tree = parse("((A,(B1,B2)B)C,D)")
    var nodes: seq[string]
    for n in traverseBF(tree.nodes[4]):  # This is the node labeled "C"
      nodes.add n.label
    doAssert nodes == @["C", "A", "B", "B1", "B2"]

  var q = [n].toDeque
  while len(q) > 0:
    let v = q.popFirst
    yield v
    for child in v.children:
      q.addLast child

iterator traverseDF*(n: Node): Node =
  ## Traverses (searchs) a tree starting at node `n` and yielding nodes using
  ## `Depth First Search`_.
  ##
  ## .. _`Depth First Search`: https://en.wikipedia.org/wiki/Depth-first_search
  ##
  ## E.g. Given the tree:
  ##
  ## .. Don't delete the empty line below (it's a \u00a0 char), or the tree
  ##    doesn't render
  ## .. code-block::
  ##
  ##   ??
  ##            ??????A1
  ##          ???A???
  ##          ??? ??????A2
  ##        ???C??? ??????B1
  ##      ????????? ???B???
  ##        ???   ??????B2
  ##        ??????????????????D
  ##
  ## Searching Depth First from node C, it yields B, B2, B1, A, A2 and A1, in
  ## that order.
  runnableExamples:
    let tree: Tree = parse("((A,(B1,B2)B)C,D)")
    var nodes: seq[string]
    for n in traverseDF(tree["C"]):
      nodes.add n.label
    doAssert nodes == @["C", "B", "B2", "B1", "A"]

  var s = @[n].toDeque
  while len(s) > 0:
    let v = s.popFirst
    yield v
    for child in v.children:
      s.addFirst child

proc parse*(s: string): Tree =
  ## Parses a string into a `Tree<newick.html#Tree>`_.
  ##
  runnableExamples:
    let tree: Tree = parse("((A,B),C)")
    doAssert len(tree) == 5  # 3 leafs + 1 branch (A,B) + 1 branch (A,B),C
    doAssert tree.nodes[0].label == "A"

  let seps = {'(', ')', ',', ':', ';', '\'', '[', ']'}
  result = Tree()
  var nodes: seq[Node] = @[Node()]

  var isComment, isLength, isQuoted: bool = false
  var label: string

  for token in tokenize(s, seps):
    if token.isSep:
      # Seps come grouped: e.g. ??,(?? or ??(((??
      # We have to consume the separators and take decissions
      var seps = token.token
      if isQuoted and seps.startsWith("''"):
        # Add double single quotes as single quotes to quoted labels.
        nodes[^1].label.add "'"
        seps = seps[2 ..^ 1]

      for sep in seps:
        if isComment and sep != ']':
          # Add anything in brackets to the comment
          nodes[^1].comment.add sep
          continue
        if isQuoted and sep != '\'':
          # If we are in a quoted label, don't process seps unless it's ??'??
          label.add seps
          continue
        case sep:
        of '(':
          nodes.add Node(parent: option(nodes[^1]))
          nodes[^2].children.add(nodes[^1])
        of ')', ';':
          result.nodes.add nodes.pop
        of ',':
          result.nodes.add nodes.pop
          nodes.add Node(parent: result.nodes[^1].parent)
          nodes[^2].children.add(nodes[^1])
        of ':':
          isLength = true
        of '\'':
          isQuoted = not isQuoted
        of '[':
          isComment = true
        of ']':
          isComment = false
        else:
          # This shouldn't happen
          discard
    else:  # This is a value (a label, a comment or a length)
      let value = collect(newSeq):
        for c in token.token:
          if not isQuoted and not isComment:
            if c == '_':
              ' '
            elif c notin Whitespace:
              c
          elif c notin Newlines:
            c

      if isComment:
        nodes[^1].comment.add join(value)
      elif isLength:
        nodes[^1].length = option(parseFloat(join(value)))
        isLength = false
      else:
        if nodes.len > 0:
          nodes[^1].label.add label & join(value)
          label = ""

  if len(nodes) == 1: result.nodes.add nodes.pop

proc parse*(strm: Stream): Tree =
  ## Parses a stream into a `Tree<newick.html#Tree>`_.
  ##
  parse(strm.readAll)
