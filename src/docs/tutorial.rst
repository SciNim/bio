=================
Bio: the tutorial
=================

:Author: Xabier Bello
:Version: |libversion|
:Nim Version: |nimversion|

.. contents::


Creating sequences
==================


.. code-block:: nim

    import bio/sequences

    let mySeq = Sequence(chain: "ACGTGTGAAC")


.. code-block:: nim

    import bio/sequences

    let myDNA = guess("ACGTGTGAAC")

    echo myDNA.complement

API: `sequences <sequences.html>`_.

Working with files
===================

Currently you can use only Fasta files. There are a number of ways of dealing
with Fasta files:

* Loading the whole file in memory. Only if the file is small.
* Iterating the whole file while doing whatever you need to each sequence (e.g.
  GC counting, length counting, sequence padding...), and close the file.
* Opening a file and selecting one or a few sequences to do something.

Loading the whole file
----------------------

Download this file_ to a directory of your choosing, and uncompress it to
follow the examples. In the same directory open a new file 'tutorial1.nim' in
your code editor and type:

.. _file: https://gitlab.com/xbello/bio/-/blob/master/src/docs/sample.fasta.gz

.. code-block::

   import bio/fasta

   let srs: seq[SequenceRecord] = load("sample.fasta")
   for sequence in srs:
     echo sequence.name, " GC-content: ",
       sequence.record.chain.count({'C', 'G'}) /
       sequence.record.chain.count({'A', 'C', 'G', 'T'})

Exit the editor and compile/run the code:

.. code-block::

    nim c -r tutorial1.nim

You should see some messages from the compiler like "Hint: system [Processing]"
and after that a lot of lines with the GC-content for each sequence:

.. code-block::

  AF517068.1 HIV-1 isolate 296-6 from Spain retrotranscriptase (pol) gene,\
  partial cds GC-content: 0.3668261562998405

Basically, 'srs' is a sequence (in other languages, a *list*) of
`SequenceRecord` objects. Each of this objects has at least a `name` and an
inner `Sequence` as `record` that has the chain as a string. You can use any
procedure you'll use with a string, such the ones in strutils_, on this chain.

.. _strutils: https://nim-lang.org/docs/strutils.html

You can also use the methods for `Sequences` we have seen above. Modify
'tutorial1.nim' to this:

.. code-block::

   import bio/fasta

   let srs: seq[SequenceRecord] = load("sample.fasta")
   for sr in srs:
     echo sr.sequence.translate

.. code-block::

   Protein: FVQIWKRKGKFQKLGLKIHIILQCLP*RKKTVQNGEN**ISENLIRKLKTSGKFN*EYHI…
   ...
   Protein: YASRQLKTHERNNTMHDLELGAVVFALKIWRHYLYGTRCTIYTDHMSLEHIFKQKKLNIR…
   Protein: DVKTAFLHGVVEEEVYVEQPPGLEDPLHPDRVLLLNKALYGHHQAPRAWYETLSTYLLSN…

Some of the sequences seems to be in frame, while others are not. Let's try a
naive method to find the longest ORFs of each sequence:

.. code-block::

  import strutils
  import bio/fasta

  let srs: seq[SequenceRecord] = load("sample.fasta")

  var currentTrans, bestTrans: Sequence
  var stops, currentStops: int

  for sr in srs:
    stops = sequenceRecord.len
    for shift in 0 .. 2:
      currentTrans = sr[shift .. ^(shift + 1)].record.translate
      currentStops = currentTrans.chain.count('*')
      if currentStops < stops:
        stops = currentStops
        bestTrans = currentTrans
    echo bestTrans, " STOPS: ", stops

In this code we create some variables (we will mutate it multiple times) in
lines 6 and 7. In line 9 we cycle every `SequenceRecord` loaded before, givin
`stops` an initial value of the maximum number of `*` possible.

Line 11 cycles through each of the three frameshifts in the direct sequence,
translating each one in line 12, counting the number of `*` found in line 13
and storing this sequence as the best of the three (so far) if the number of
stops is lower than before in lines 14-16.

In the line 17 we echo the `Protein Sequence` and the number of stops.

    If you are having problems with Nim's notation for sequence slicing, go
    take a tour at the `slice Tutorial`_.
    Always leave spaces between the indexes and the inner operator (`..`).

.. slice tutorial_: https://nim-lang.org/docs/tut1.html#advanced-types-slices>

Sequence modifying
------------------

`Bio` allows you to modify the sequences in place. Use this power with care, as
we don't have undo buttons here. Let us try with another example to pad the
sample file with `-` until they are all equal length, as some software doesn't
like files with uneven sequences:

.. code-block::

  import strutils
  import bio/fasta

  let srs: seq[SequenceRecord] = load("sample.fasta")

  var maxLen: int
  for sr in srs:
    maxLen = max(len(sr), maxLen)

  for sr in srs:
    sr.record.chain.add repeat('-', maxLen - len(sr))

  srs.dumpTo("sample_output.fasta")

Wow, that was easy. In a first loop we find the longest sequence of the batch
(lines 6-8). Then we modify *in place* every sequence adding as many `-` as
they need to reach `maxLen` value.

In the last line we save, or dump, the data into a new file. If you copy-pasted
the name of the file from the line 4 to the line 13, your original fasta will
be overwritten without warning.

    `dumpTo` is called so because I'm copying the `Python JSON`_ module naming.
    The `To` reminds me that "dump sequences to file", because I was prone to
    read "file receives dumped sequences". The Nim marshal_ module uses
    `store` instead of `dump` and maybe it's a better naming, while it seems to
    reserve `dump` to `thing in memory gets dumped into screen`.

.. _Python Json: https://docs.python.org/3/library/json.html
.. _marshal: https://nim-lang.org/docs/marshal.html

* `fasta <fasta.html>`_ is a set of utilities to work with Fasta files.
