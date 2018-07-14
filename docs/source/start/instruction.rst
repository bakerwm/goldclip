Getting Started
================

The idea of developing ``sub-commands`` in Python was from an package end2you_ on Github.

Recently, I was working on a project to develop a pipline for bioinformatics analysis, which is consist of several main steps, eg: trimming, mapping, counting, DE analysis. For actual analysis, we might need to invoke the pipeline from any of the step, or we just do mapping.

.. _end2you: https://github.com/end2you/end2you.git

Installation
-------------

Download this repo from github and enter the directory:

::

    $ git clone https://github.com/bakerwm/demo.git
    $ cd demo
    $ python setup.py install 

You can test the installation by ``demo -h``. It is ok if you can see the following message:

::

    $ demo -h
    usage: demo [-h] {map,count,diff} ...

    Group of sub-commands.

    positional arguments:
      {map,count,diff}  sub-commands.
        map             mapping reads to reference.
        count           count reads on features.
        diff            differentially expression analysis.

    optional arguments:
      -h, --help        show this help message and exit


Overview
---------

Structure of the directory

::

    .
    ├── demo
    │   ├── demo.py
    │   ├── __init__.py
    │   └── programs
    │       ├── count.py
    │       ├── diff.py
    │       ├── __init__.py
    │       ├── map.py
    │       └── parser.py
    ├── docs
    │   ├── build
    │   ├── Makefile
    │   └── source
    │       ├── conf.py
    │       ├── index.rst
    │       ├── start
    │       │   └── instruction.rst
    │       ├── _static
    │       └── _templates
    ├── LICENSE
    ├── README.md
    └── setup.py

Details
--------

The main components of the packages are: ``setup.py``, ``demo/`` and also instruction files ``LICENSE`` and ``README.md``.

+ ``demo/demo.py``

It is the main entery for the package, parsing arguments, choosing which sub-command to invoke.

+ ``demo/programs/``

Saving the sub-commands: ``map.py``, ``count.py``, ``diff.py`` and ``parser.py``

+ ``__init__.py``

In order to make sub-commands available, add something in the ``__init__.py`` works.

::

    $ cat demo/__init__.py
    from .programs import *

    $ cat demo/programs/__init__.py
    from .map import *
    from .count import *
    from .diff import *
    from .parser import *


FAQ
----

+ 1. How to use other languages besides Python?

+ 2. How to use extra data within the package?














