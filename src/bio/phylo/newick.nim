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
##       echo node.label, ">", nodel.length
##

import options, streams, strutils, sugar


# TODO: types should go to some phylo base file
type
  Node* = ref object
    label*, comment*: string
    length*: Option[float]
    children*: seq[Node]
    parent*: Node

  Tree* = ref object
    nodes*: seq[Node]

proc len*(t: Tree): Natural =
  len(t.nodes)

proc parse*(s: string): Tree =
  ## Parses a string into a `Tree<newick.html#Tree>`_.
  ##
  runnableExamples:
    let tree: Tree = parse("((A,B),C)")
    echo len(tree)
    doAssert len(tree) == 5  # 3 leafs + 1 branch (A,B) + 1 branch (A,B),C
    doAssert tree.nodes[0].label == "A"

  let seps = {'(', ')', ',', ':', ';', '\'', '[', ']'}
  result = Tree()
  var nodes: seq[Node] = @[Node()]

  var isComment, isLength, isQuoted: bool = false
  var label: string

  for token in tokenize(s, seps):
    if token.isSep:
      # Seps come grouped: e.g. «,(» or «(((»
      # We have to consume the separators and take decissions
      var seps = token.token
      if isQuoted and seps.startsWith("''"):
        # Add double single quotes as single quotes to quoted labels.
        nodes[^1].label.add "'"
        seps = seps[2..^1]

      for sep in seps:
        if isComment and sep != ']':
          # Add anything in brackets to the comment
          nodes[^1].comment.add sep
          continue
        if isQuoted and sep != '\'':
          # If we are in a quoted label, don't process seps unless it's «'»
          label.add seps
          continue
        case sep:
        of '(':
          nodes.add Node(parent: nodes[^1])
          nodes[^2].children.add(nodes[^1])
        of ')':
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
        of ';':
          result.nodes.add nodes.pop
        else:
          # This shouldn't happen
          discard
    else:  # This is either a label, a comment or a length
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
