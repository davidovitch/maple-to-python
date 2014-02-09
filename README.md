maple-to-python
===============

Introduction
------------

This is a very early and naive prototype on an incomplete and unfinished 
conversion scheme from a Maple [0] worksheet to either an IPython notebook or a
Python/SymPy script.

As a first step, the Maple commands are extracted from the worksheet (*.mw) XML 
file, and converted to text using mw2txt.py [1](http://blog.tremily.us/posts/Maple/).
The second step takes the text output of the Maple worksheet and uses pyparsing [2]
to parse the commands and convert them to SymPy [3] compatible commands.

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

TODO:

* Merge mw2txt.py and mw2py.py into a single library

* not all Maple function calls are translated correctly to SymPy constructs.
Either write wrapper functions with the same Maple syntax and includ them in the
output, or create additional parsing rules to translate them to native SymPy
calls.

* Write a lot of unittests to make sure the parsing is done correctly for as
meany corner cases as possible

* Save as IPython notebook or Python script

[0] [Maple](http://www.maplesoft.com/products/maple/)
[1] [mw2txt.py](http://blog.tremily.us/posts/Maple/)
[2] [pyparsing](http://pyparsing.wikispaces.com/)
[3] [SymPy](http://sympy.org/en/index.html)

Maple is a registered trademark of Maplesoft.

