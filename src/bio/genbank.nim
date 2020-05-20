## :Author: |author|
## :Version: |libversion|
import strutils
import strtabs

import sequences
export sequences


const FIELD = 12

proc parseFeatures(fileIn: File): seq[Feature] =
  ## Private procedure
  ##
  ## Passing a file with the reading pointer at the "FEATURES" line, this proc
  ## parses all the Features found and return a seq[Feature] to the caller.
  ##
  ## The file pointer after the proc is at the next Field (usually ORIGIN).
  var line = fileIn.readline
  var lineData: seq[string]
  var thisFeature: Feature
  var currentQualifier: string
  var qualifiers: StringTableRef

  while line[0] == ' ':
    if line[5] != ' ':  # This is a new key line
      if thisFeature != nil:  # Save the previous feature, init another one
        thisFeature.qualifiers = qualifiers
        result.add thisFeature
      qualifiers = newStringTable(modeCaseSensitive)
      thisFeature = Feature(key: line[5 .. 20].strip(),
                            location: line[21 .. ^1].strip(),
                            qualifiers: qualifiers)
    else:
      lineData = line[22 .. ^1].split("=")
      if line[21] == '/':  # This is a new qualifier
        currentQualifier = lineData[0]
        if currentQualifier in qualifiers:
          # Some qualifiers can appear multiple times (usually db_xref)
          qualifiers[currentQualifier].add "," & lineData[1].replace("\"", "")
        else:
          qualifiers[currentQualifier] = lineData[1].replace("\"", "")
      else:  # This is a continuation line from the previous qualifier
        qualifiers[currentQualifier].add lineData[0].replace("\"", "")
    line = fileIn.readline

  # This is the last data that remains hanging in variables
  qualifiers[currentQualifier].add lineData[0].replace("\"", "")
  result.add thisFeature

proc load*(fName: string, kind: string="gb"): seq[SequenceRecord] =
  ## Load a `seq` of `SequenceRecords<sequences.html#SequenceRecord>`_ from a
  ## filename.
  ##
  ## *Careful*: this loads **all** the file at once. If you need an interator
  ## over the file, use the `sequences iterator<#sequences.i,string,string>`_
  ##
  ## .. code-block::
  ##
  ##    import bio/genbank
  ##
  ##    let mySeqs = load("path/to/file.gb")
  let fileIn: File = open(fName)
  defer: fileIn.close

  var line, name, sequence: string
  var features: seq[Feature]

  while not fileIn.endOfFile:
    line = fileIn.readLine
    if line.startsWith("LOCUS"):
      discard ## TODO Get the sequence type and the locus name
    if line.startsWith("ACCESSION"):
      discard ## TODO Get the accession id
    if line.startsWith("DEFINITION"):
      name = line[FIELD .. ^1]
      line = fileIn.readLine
      while line[0] == ' ':
        name.add line[FIELD - 1 .. ^1]  # (- 1) borrows a space to join lines
        line = fileIn.readLine
    if line.startsWith("FEATURES"):
      features = parseFeatures(fileIn)
    if line.startsWith("ORIGIN"):
      line = fileIn.readLine
      while line[0] == ' ':
        sequence.add line[10 .. ^1].replace(" ", "")
        line = fileIn.readLine
    if line.startsWith("//"):
      result.add SequenceRecord(name: name,
                                record: guess(sequence.toUpperAscii),
                                features: features)
