What is Bio
===========

<img src="src/media/logo.svg" width="200" />
A library to work with biological sequences. Using Nim.

How to install
==============

In you have Nim and Nimble installed, type (version is optional):

    $ nimble install https://github.com/SciNim/bio 0.2.7

Click here if you need to install [Nim](https://nim-lang.org/).

How to get/read the docs
========================

Already built docs are at https://xbello.gitlab.io/bio/

You can also build your own docs.  With the nimble tool, type at the root of
the repo:

    $ nimble repodocs

The script creates a "public" directory, scans the source code for docs,
runs the embedded doc-tests and finally build the documentation as a webpage.
Open the file "public/index.html" in any browser.

The "public" directory can be moved around to any place of your convenience.

How to work with bio
====================

This library doesn't provide any functionality as executables out of the box.
Refer to the Tutorial and Cookbook in the Documentation to get a glimpse of
what you can do with it.
