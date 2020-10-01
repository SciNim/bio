==============
Bio: a library
==============

:Author: Xabier Bello
:Version: |libversion|
:Nim Version: |nimversion|

.. contents::


Install
=======

To install the library the easiest way is to use:

.. code-block::

    $ nimble install bio

As this is a library, you should include `bio` in your .nimble file. I always
prefer to pin the versions, as the API may change between releases and break
your code. In your .nimble you will find the following lines:

.. code-block:: nim

    # Dependencies

    requires "nim >= 1.2.0"

Just add a new requirement or expand the existing one:

.. code-block:: nim

    requires "nim >= 1.2.0"
    requires "bio 0.2.0"

.. code-block:: nim

    requires "nim >= 1.2.0, bio 0.2.0"

..

    **Q:** Is there any Docker image or equivalent? No.

Build the docs
--------------

Once you installed the library, the docs are ready to be created locally.
Nimble install the packages at ~/.nimble/pkgs/bio-|libversion| (UNIX) or at
[TBD: windows path]. Go there and execute:

.. code-block::

    $ nimble docs

Doctests will be run and if everything goes smoothly, you will get a local copy
of the docs at `htmldocs/`, in sync with your version, that you can browse
opening `htmldocs/index.html`. This path is self-contained, i.e. you can move
it wherever you like and it will work. Even if the library updates and the docs
changes, your local copy will be intact.

If you prefer the latests online version, go to [TBD: online docs].

Usage
=====

Once `installed<#install>`_, import the parts of the library you will need
for your code:


.. code-block:: nim

    import bio/sequences
    import bio/fasta

Some modules already offer other modules (via `export`) so you don't need to
re-import them. E.g. `import bio/fasta` exports `bio/sequences`, so you can
use Sequence_, SequenceRecord_ and others without explicitly importing the
module `bio/sequence`.

Speed considerations
--------------------

Chances are you are using this library instead something like
`Biopython <http://biopython.org/>`_ because you expect some speed added to
your code for free, because "Nim is faster than Python". This will **not** be
the case. Biopython is fast enough, and in fact you have to pay attention to
your Nim-fu if you don't want to end with slower code.

An example, scanning sequences in a 500 Mb Fasta file with 5000 sequences:

..example-code::

  .. code-block:: nim

      # Nim code

      let seqs = @["A", "B", "C", "D", "E"]

      for sequence in sequences("500Mb_Database.fasta"):
        if sequence.name in seqs:
          echo sequence

  .. code-block::

      # Python code using Biopython

      seqs = ["A", "B", "C", "D", "E"]

      for record in SeqIO.parse("500Mb_Database.fasta", "fasta"):
          if record.id in seqs:
              print(record.id)


============================================   ===========
  Command                                         Time
============================================   ===========
`nim c program.nim`                               20.226 s
`nim c -d:release program.nim`                     2.466 s
`python program.py`                                1.634 s
`nim c -d:danger program.nim`                      1.580 s
`nim c -d:danger --passC:-flto program.nim`        1.074 s
============================================   ===========

The same simple code could give you a 10x slower or a 50% faster program than
Python, just by changing the compilation options.

Much more powerful reasons to use Nim + bio would be the distribution of
binaries (you should be safe compiling in your computer and then moving the
program around as a single file), the FFI_ (easy access to libraries like hts_)
or the ease to do multithread_.


.. _FFI: https://nim-lang.org/docs/manual.html#foreign-function-interface
.. _multithread: https://nim-lang.org/docs/manual.html#threads
.. _hts: http://www.htslib.org/


Tutorial / Recipes
==================

`Tutorial <tutorial.html>`_ is a gently introduction to the API. Once you get
a grasp on the basics, `Recipes <recipes.html>`_ includes some snippets that
mix together elements of the API to get some task done.

API docs
========

The API docs include all the `Objects`, `procs` and companions documented in
isolation. I tried to include code samples where I could, but I find quite
difficult to learn something going straight to the API docs. It should be your
main reference once you get a bit familiar with the library through the
`Tutorial <tutorial.html>`_.

Sequences
---------

* `sequences <sequences.html>`_ explains Sequences, the core of the library.

Operations with files
---------------------

* `fasta <fasta.html>`_ is a set of utilities to work with Fasta files.

NCBI Entrez
-----------

* `entrez <entrez.html>`_ allows to interact with the Entrez tools.

.. _Sequence: sequences.html#Sequence
.. _SequenceRecord: sequences.html#SequenceRecord
