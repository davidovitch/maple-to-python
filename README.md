maple-to-python
===============

Introduction
------------

This is a very early and naive prototype on an incomplete and unfinished 
conversion scheme from a [Maple](http://www.maplesoft.com/products/maple/) 
worksheet to either an IPython notebook or a Python/SymPy script.

As a first step, the Maple commands are extracted from the worksheet (*.mw) XML 
file, and converted to text using [mw2txt.py](http://blog.tremily.us/posts/Maple/).
The second step takes the text output of the Maple worksheet and uses 
[pyparsing](http://pyparsing.wikispaces.com/) to parse the commands and convert
them to [SymPy](http://sympy.org/en/index.html) compatible commands.

mw2txt.py
---------

Originally written by W. Trevor King in order to display Maple worksheets in a
shell, optionally outputted in color. For example:

```
$ mw2txt.py example.mw 
Hi there
> restart;
> interface(prettyprint=0):
> 1;# one  + plus 2 two ;
1
> 3 + 4;  bold
7
Equation
```

Roadmap
-------

* Merge mw2txt.py and mw2py.py into a single library

* not all Maple function calls are translated correctly to SymPy constructs.
Either write wrapper functions with the same Maple syntax and include them in the
output, or create additional parsing rules to translate them to native SymPy
calls.

* The Maple syntax can be challenging to interpreted: the square brackets
are used either to create variable names with subscripts, or when a vector, it is 
used to indicate index of that vector. A more robust approach needs to be
implemented to correctly distinguish vectors from variables with subscripts.

* Write a lot more unittests to make sure the parsing is done correctly for as
many corner cases as possible

* Save as [IPython notebook](http://ipython.org/notebook.html) or Python script

Maple is a registered trademark of Maplesoft.

