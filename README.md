Some of the classes I TA use [Maple][].  (Caveat: I prefer [[Python]],
as a more general language.  Use [[SymPy]] or [Sage][] if you need
symbolic processing.)  Anyhow, I get Maple worksheets to grade.
[[SSH]]ing into the department computer lab to fire up `xmaple` is a
pain, so I wrote [[mw2txt.py]] to extract the Maple commands from the
worksheet.  It benefits from the fact that worksheets are fairly clean
XML.  Graphs and equations are more difficult, since they have
complicated layout and are stored as encoded blobs.  Other than that,
things work pretty well.  Here's the output from my [[example.mw]]
example worksheet, picking out the math-mode sections (in red) and
unprocessed blocks (in yellow) from the comments (uncolored).

<pre><code>$ mw2txt.py --color example.mw 
Hi there
<span style="color: red">&gt; restart;
&gt; interface(prettyprint=0):
&gt; 1;#</span> one  <span style="color: red">+</span> plus <span style="color: red">2</span> two <span style="color: red">;</span>
1
<span style="color: red">&gt; 3 + 4;</span>  bold
7
<span style="color: yellow">Equation</span>
</code></pre>

[Maple]: http://www.maplesoft.com/products/maple/
[Sage]: http://www.sagemath.org/

[[!tag tags/code]]
[[!tag tags/python]]
