#!/usr/bin/env python
# 
# Copyright (C) 2011-2012 W. Trevor King <wking@drexel.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

"""View Maple worksheets (.mw) from the command line (without X).

./mw2txt.py -c input.mw | less
./mw2txt.py -m input.mw | /opt/maple15/bin/maple | less
"""

import sys as _sys

import lxml.etree as _lxml_etree

try:
    from pygments.console import colorize as _colorize
except ImportError, e:
    _sys.stderr.write(str(e) + '\n')
    def _write_color(string, color=None, stream=None):
        if stream is None:
            stream = _sys.stdout
        stream.write(string)
else:
    def _write_color(string, color=None, stream=None):
        if color is None:
            color = 'reset'
        if stream is None:
            stream = _sys.stdout
        stream.write(_colorize(color_key=color, text=string))


__version__ = '0.1'


class Writer (object):
    def __init__(self, color=None, stream=None, use_color=False):
        self.color = color
        self.stream = stream
        self.use_color = use_color
        self.last_char = None

    def __call__(self, text, color=None):
        if not self.use_color:
            color = None
        elif color is None:
            color = self.color
        if text == '\n' and self.last_char == '\n':
            return  # don't add lots of blank lines
        _write_color(string=text, color=color, stream=self.stream)
        self.last_char = text[-1]


def mw2txt(path, writer, filter_math=False):
    xml = _lxml_etree.parse(path)
    pruned_iteration(
        root=xml.getroot(),
        match=lambda node: node.tag == 'Text-field',
        match_action=lambda node: top_text_node2txt(
            node=node, writer=writer, filter_math=filter_math),
        match_tail=lambda node:writer(text='\n'))

def top_text_node2txt(node, writer, filter_math=False):
    if filter_math:
        match_action = None
    else:
        match_action = lambda node: other_in_text_node2txt(
            node=node, writer=writer)
    pruned_iteration(
        root=node,
        match=lambda node: node.tag not in ['Text-field', 'Font', 'Hyperlink'],
        match_action=match_action,
        other_action=lambda node: text_node2txt(
            node=node, writer=writer, filter_math=filter_math),
        match_tail=lambda node:tail_node2txt(
            node=node, writer=writer, filter_math=filter_math),
        other_tail=lambda node:tail_node2txt(
            node=node, writer=writer, filter_math=filter_math))

def other_in_text_node2txt(node, writer):
    if node.tag in ['Drawing-Root']:
        # ignore missing content
        pass
    elif node.tag in ['Equation', 'Image', 'Plot']:
        # warn about missing content
        writer(text=node.tag, color='yellow')
    else:
        # warn about wierd tag
        writer(text=node.tag, color='magenta')

def text_node2txt(node, writer, filter_math=False):
    if node.tag not in ['Text-field', 'Font', 'Hyperlink'] and not filter_math:
        # warn about wierd tag
        writer(text=node.tag, color='magenta')
    write_text(
        node=node, text=node.text, writer=writer, filter_math=filter_math)

def tail_node2txt(node, writer, filter_math=False):
    if node.tag != 'Text-field':
        write_text(
            node=node.getparent(), text=node.tail, writer=writer,
            filter_math=filter_math)

def write_text(node, text, writer, filter_math=False):
    if not text:
        return
    style = node_style(node)
    if filter_math:
        if style == 'Maple Input':
            writer(text=text)
        return
    prompt = node.get('prompt', None)
    if prompt:
        t = '\n'.join(prompt+line for line in text.splitlines())
        if text.endswith('\n'):
            t += '\n'  # '\n'.join('a\nb\n'.splitlines()) == 'a\nb'
        if writer.last_char not in [None, '\n']:
            t = t[len(prompt):]  # no initial prompt
        text = t
    if style == 'Maple Input':
        color = 'red'
    else:
        color = None
    writer(text=text, color=color)

def node_style(node):
    p = node
    while p is not None:
        style = p.get('style', None)
        if style:
            return style
        p = p.getparent()
    return None

def pruned_iteration(root, match, match_action=None, match_tail=None,
                     other_action=None, other_tail=None, debug=False):
    if debug:
        _write_color('start pruned iteration from %s\n' % root, color='blue')
    line = [None]
    stack = [root]
    while len(stack) > 0:
        node = stack.pop(0)
        p = node.getparent()
        while line[-1] != p:
            n = line.pop()
            if n is None:
                break
            _pruned_iteration_handle_tail(
                node=n, match=match, match_tail=match_tail,
                other_tail=other_tail, debug=debug)
        line.append(node)
        if debug:
            color_node(node, color='cyan')
        if match(node):
            if match_action:
                match_action(node)
        else:
            if other_action:
                other_action(node)
            stack = list(node.getchildren()) + stack
    while len(line) > 0:
        n = line.pop()
        if n is None:
            break
        _pruned_iteration_handle_tail(
            node=n, match=match, match_tail=match_tail, other_tail=other_tail,
            debug=debug)
    if debug:
        _write_color('end pruned iteration from %s\n' % root, color='blue')

def _pruned_iteration_handle_tail(node, match, match_tail, other_tail,
                                  debug=False):
    if debug:
        color_node(node, color='magenta', tail=True)
    if match(node):
        if match_tail:
            match_tail(node)
    else:
        if other_tail:
            other_tail(node)

def node_depth(node):
    depth = 0
    p = node.getparent()
    while p is not None:
        depth += 1
        p = p.getparent()
    return depth

def color_node(node, color=None, tail=False):
    depth = node_depth(node)
    string = ' '*depth + node.tag
    if tail:
        string += ' tail'
    _write_color(string + '\n', color)


if __name__ == '__main__':
    from optparse import OptionParser as _OptionParser

    # don't wrap epilog paragraphs
    class OptionParser (_OptionParser):
        def format_epilog(self, formatter):
            return self.epilog

    parser = OptionParser(
        usage='%prog [options] input.mw', epilog='\n'+__doc__)
    parser.add_option(
        '-c', '--color', dest='color', action='store_true',
        help='Use ANSI escape sequences to color output')
    parser.add_option(
        '-m', '--maple', dest='maple', action='store_true',
        help='output text suitable for piping into `maple`')

    options,args = parser.parse_args()
    path = args[0]

    filter_math = options.maple
    writer = Writer(use_color=options.color)
    if options.maple:
        if options.color:
            raise ValueError("maple doesn't understand ANSI color")
    mw2txt(path=path, writer=writer, filter_math=filter_math)
