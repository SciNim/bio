## :Author: |author|
## :Version: |libversion|
import strformat
import strutils
import tables

import sequences
export sequences


type
  FileType* = enum
    ## The filetype, to allow potential future filetypes handling.
    ftFasta = "fasta"

  Index* = ref object of RootObj
    ## The index is a `table` with the name of the sequence as a string, except
    ## the initial '>', and its position in the file. The filename which the
    ## Index refers to is stored in `source`.
    source*: string
    table*: TableRef[string, int64]

## Indexing vs iterating
## =====================
##
## Given a Fasta file, 1.3 Gb in size, with ~40000 sequences and 30 kb per
## sequence, do the following operations compiled with `-d:danger` and timing
## with `/usr/bin/time -v`:
##
## ====================================   =========== ===============
##  Operation                                 Time        Memory
## ====================================   =========== ===============
##  newIndex("yourFasta.fas")                 1.90 s     ~6000 kbytes
##  Index + search 1 known sequences...
##  ... at beginning of the file              2.04 s     ~6000 kbytes
##  ... at the middle                         2.10 s        "
##  ... at the end                            2.07 s        "
##  Index + search 100 random sequences       2.08 s     ~6500 kbytes
##  Index + search all sequences              4.80 s     ~6500 kbytes
##
##  Load pre-saved Index...
##  ... and search 1 known sequence           0.02 s    ~15000 kbytes
##  ... and search 100 random sequences       0.03 s        "
##
##  Iterate through all sequences             2.45 s     ~4200 kbytes
##  ====================================   ===========  =============
##
##
## Random sequence extraction
## ==========================
##
## You can get partial sequences using the `[]` procedures in this module. But
## if you plan to use a lot of querying on huge chromosomes, maybe you should
## try `hts-nim<https://github.com/brentp/hts-nim>`_. The feature is a bit
## hidden but here you have a sample (needs an `index of the fasta
## <https://www.htslib.org/doc/samtools-faidx.html>`_ file):
##
##
## .. code-block::
##
##   import hts
##
##   var faidx: Fai
##
##   doAssert f.open("UCSC.hg19.fasta")
##   echo f.get("chr1:100000-100100")
##
##   # Or
##
##   echo f.get("chr1", 100000, 100100)
##
## That is about 100x faster and uses 100x less memory than the equivalent
## `bio` code.
##
iterator sequences*(data: seq[string], kind: FileType=ftFasta):
  SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## sequence of strings, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## I found this useful only to embed a sequence in a program.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   const data = slurp("path/to/file.fas").splitLines
  ##
  ##   for sequence in sequences(data):
  ##     doAssert(sequence of SequenceRecord)
  # XXX This can be probably done unified with the other sequences(file)
  #     through Streams
  let last: int = data.len - 1
  var name, sequence: string
  for i, line in data.pairs:
    if (line.startsWith('>') and sequence.len > 0) or i == last:
      if i == last:
        sequence.add line
      yield SequenceRecord(name: name, record: guess(sequence.toUpperAscii))

    if line.startsWith('>'):
      sequence = ""
      name = line[1 .. ^1]
    else:
      sequence.add line

iterator sequences*(fName: string, kind: FileType=ftFasta):
  SequenceRecord {.inline.} =
  ## Iterate through all the `Sequences<sequences.html#Sequence>`_ in a given
  ## filename, yielding `SequenceRecords<sequences.html#SequenceRecord>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   for sequence in sequences("path/to/file.fas"):
  ##     doAssert(sequence of SequenceRecord)

  let fileIn: File = open(fName)
  defer: fileIn.close

  var name, sequence: string
  for line in fileIn.lines:
    if (line.startsWith('>') and sequence.len > 0) or fileIn.endOfFile:
      if fileIn.endOfFile:
        sequence.add line
      yield SequenceRecord(name: name, record: guess(sequence.toUpperAscii))

    if line.startsWith('>'):
      sequence = ""
      name = line[1 .. ^1]
    else:
      sequence.add line

proc newIndex*(fName: string): Index =
  ## Build an `Index<#Index>`_ for a given fasta file.
  ##
  ## You should build an index if you plan to access a few sequences by name
  ## randomly and repeatedly. If you'll need all the sequences and specially
  ## if you only need it once, use `sequences
  ## iterator<#sequences.i,string,string>`_
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##
  ## Once created, saving and loading of Indexes is covered in a
  ## `recipe<recipes.html#building-storing-and-retrieving-an-index>`_.
  ##
  var table = newTable[string, int64]()

  let fileIn: File = open(fName)
  defer: fileIn.close

  const newLine = 1

  for line in fileIn.lines:
    if line.startsWith('>'):
      table[line[1 .. ^1]] = fileIn.getFilePos - len(line) - newLine

  Index(source: fName, table: table)

proc `[]`*(index: Index, fileIn: File, name: string): SequenceRecord {.inline.} =
  ## Returns a `SequenceRecord<sequences.html#SequenceRecord>`_ from an
  ## `Index<#Index>`_ using an already opened file handler.
  ##
  ## Ignores the file path provided with `Index.source` field, and doesn't
  ## check the validity of `fileIn` handler (e.g. if it's in
  ## `fmRead<https://nim-lang.org/docs/io.html#FileMode>`_ mode).
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##
  ##   block:
  ##     let fastaFile = open("path/to/file.fas")
  ##     defer: fastaFile.close
  ##     echo index[fastaFile, "My Record"]
  ##     echo index[fastaFile, "Other Record"]

  fileIn.setFilePos(index.table[name])
  let name = fileIn.readLine[1 .. ^1]
  var sequence: string
  for line in fileIn.lines:
    if line.startsWith('>') or fileIn.endOfFile:
      if fileIn.endOfFile:
        sequence.add line
      return SequenceRecord(name: name,
                            record: guess(sequence.toUpperAscii))
    sequence.add line

proc `[]`*(index: Index, name: string): SequenceRecord {.inline.} =
  ## Returns a `SequenceRecords<sequences.html#SequenceRecord>`_ from an
  ## `Index<#Index>`_.
  ##
  ## Doesn't check if the file provided at `Index.source` exists, and
  ## open/closes it on each access.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##   echo index["My Record"]

  let fileIn: File = open(index.source)
  defer: fileIn.close

  index[fileIn, name]

proc `$`*(index: Index): string {.inline.} =
  ## Returns a brief description of the `Index<#Index>`_.
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##   echo index
  ##
  ## .. code-block::
  ##
  ##   > Index for path/to/file.fas, length: 232

  &"Index for {index.source}, length: {len(index.table)}"

iterator `items`*(index: Index): string =
  ## Iterate through all the `Index<#Index>`_ items, yielding the names of
  ## the sequences in the index.
  ##
  ## If you need the actual `SequenceRecords<sequences.html#SequenceRecord>`_,
  ## use the `[]` getters.
  ## If you need the positions of each SequenceRecord in file, access the
  ## underlying `Table` directly (you have to `import tables` to do it).
  ##
  ## .. code-block::
  ##
  ##   import bio/fasta
  ##
  ##   let index: Index = newIndex("path/to/file.fas")
  ##
  ##   for sequenceName in index:
  ##     echo sequenceName
  ##     echo index[sequenceName]

  for k in index.table.keys:
    yield k

proc load*(fName: string, kind: FileType=ftFasta): seq[SequenceRecord] =
  ## Load a `seq` of `SequenceRecords<sequences.html#SequenceRecord>`_ from a
  ## filename.
  ##
  ## *Careful*: this loads **all** the file at once. If you need an interator
  ## over the file, use the `sequences iterator<#sequences.i,string,string>`_
  ##
  ## .. code-block::
  ##
  ##    import bio/fasta
  ##
  ##    let mySeqs = load("path/to/file.fas")

  for sequence in sequences(fName):
    result.add sequence

proc dumpTo*(record: SequenceRecord, fHandler: File, kind: FileType=ftFasta) =
  ## Write a `SequenceRecords<sequences.html#SequenceRecord>`_ to `fHandler`,
  ## wrapping the sequence by 60 positions.
  ## The name of the SequenceRecord remains untouched.
  ##
  ## TBD: `kind` is a string as in "fasta", to support different formats.
  ## Right now only FASTA files are supported.
  ##
  ## To write a FASTA file `myOutput.fasta` with the contents:
  ##
  ## .. code-block::
  ##   `import bio/fasta`
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)
  ##
  ##   block:
  ##     let fastaOut = open("myOutput.fasta", fmWrite)
  ##     defer: fastaOut.close
  ##     myRec.dumpTo(fastaOut)
  ##
  const wrapSize: int = 60
  fHandler.write(">", record.name)

  for i, base in record.record.chain.pairs:
    if i mod wrapSize == 0:
      fHandler.write('\n')
    fHandler.write(base)
  fHandler.write('\n')
  fHandler.flushFile

proc dumpTo*(records: seq[SequenceRecord], fHandler: File, kind: FileType=ftFasta) =
  ## A shortcut to avoid the explicit cycle to write a `seq` of
  ## `SequenceRecords<sequences.html#SequenceRecord>`_
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeqA = newDna("TGCACCCCA")
  ##   let mySeqB = newDna("GTGAGAGTG")
  ##   let myRecA = SequenceRecord(name: "My DNA sequence", record: mySeqA)
  ##   let myRecB = SequenceRecord(name: "My DNA sequence", record: mySeqB)
  ##
  ##   block:
  ##     let fastaOut = open("myOutput.fasta", fmWrite)
  ##     defer: fastaOut.close
  ##     @[myRecA, myRecB].dumpTo(fastaOut)
  ##
  for sr in records:
    sr.dumpTo(fHandler, kind)

proc dumpTo*(record: SequenceRecord, fName: string, kind: FileType=ftFasta) =
  ## Same as `write-through-handler proc<#dumpTo,SequenceRecord,File,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeq = newDna("TGCACCCCA")
  ##   let myRec = SequenceRecord(name: "My DNA sequence", record: mySeq)
  ##
  ##   myRec.dumpTo("myOutput.fasta")
  ##
  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  record.dumpTo(fHandler, kind)

proc dumpTo*(records: seq[SequenceRecord], fName: string, kind: FileType=ftFasta) =
  ## Same as `write-through-handler proc<#dumpTo,seq[SequenceRecord],File,string>`_
  ## but you only need to point out the name of the file.
  ##
  ## If the file exists, it will be silently overwritten.
  ##
  ## .. code-block::
  ##   import bio/fasta
  ##
  ##   let mySeqA = newDna("TGCACCCCA")
  ##   let mySeqB = newDna("GTGAGAGTG")
  ##   let myRecA = SequenceRecord(name: "My DNA sequence", record: mySeqA)
  ##   let myRecB = SequenceRecord(name: "My DNA sequence", record: mySeqB)
  ##
  ##   myRecs = @[myRecA, myRecB]
  ##   myRecs.dumpTo("myOutput.fasta")
  ##
  let fHandler: File = open(fName, fmWrite)
  defer: fHandler.close()

  for record in records:
    record.dumpTo(fHandler, kind)
