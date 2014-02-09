# -*- coding: utf-8 -*-
"""

"""

__author__ = 'David Verelst'
__license__ = 'GPLv3'


from pyparsing import nums, Word, Literal
from pyparsing import Forward, Group, Optional
from pyparsing import alphas, ZeroOrMore, OneOrMore, ParseException
import unittest

from pprint import pprint
from collections import Iterable

def flatten(l):
    """
    source:
    http://stackoverflow.com/questions/2158395/
    flatten-an-irregular-list-of-lists-in-python/2158532#2158532
    """
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el

var_stack = {}

def flatten_to_string(iterable):
    """
    Flatten the multi-depth iterable and return as a single string
    """
    parsed_line = ''
    for k in flatten(iterable):
        parsed_line += k
    
    return parsed_line

def add_var( strg, loc, toks ):
#    var_stack.append( [f for f in flatten(toks)] )
#    var_stack.append( flatten_to_string(toks) )
    var_stack[flatten_to_string(toks)] = strg

class MapleGrammer:
    
    def __init__(self):
        
        self.bnf = self.BNF()
    
    def BNF(self):
        
        lsbr_ = Literal("[").suppress()
        rsbr_ = Literal("]").suppress()
        lsbr = Literal("[")
        rsbr = Literal("]")
        lpar  = Literal( "(" )
        rpar  = Literal( ")" )
        varname = Word(alphas + nums + "_-+.")
        equal = Literal('=')
        comma = Optional(Literal(','))
        op_opt = Optional(Word('+-*/'))

        # nested structures
        variable = Forward()
        terms_par = Forward()
        mixed_terms = Forward()
        command = Forward()
        matrix = Forward()

        # variables have nested structures with varname[subvarname[varname]]
        subname = ZeroOrMore(Group(lsbr_ + variable + rsbr_))
        variable << Group( (varname|nums) + subname )
        # keep track of all the variables added
        variable.setParseAction( add_var )
        
        # same with commands, matrices: they can be nested.
        # mind the order: first command, than variable. The variable definition
        # will be satisfied just before reaching the opening bracket and that
        # will result in a false positive for the variable
        arguments = OneOrMore( Group( (matrix|command|variable|mixed_terms) + comma) )
        command << Group( (varname + Group(lpar + arguments + rpar)) )
        
        # matrices are also nested structures
        matrix << Group(lsbr + arguments + rsbr)
        
        # terms in parenthesis can hold other nested structures
        terms = OneOrMore( Group( (variable|terms_par) + op_opt) )
        terms_par << OneOrMore( Group(lpar + terms + rpar + op_opt) )
        mixed_terms << OneOrMore( (terms_par|terms) + op_opt )
        
#        # works only for one depth
#        terms = OneOrMore( Group( (variable) + op_opt) )
#        terms_par = OneOrMore( Group(lpar + terms + rpar + op_opt) )
#        group_terms = OneOrMore( (terms_par|terms) + op_opt )
        
        # an expression of the form
        exp = Group( command | mixed_terms | variable)
        # and finally we can assemble an equation
        equation = (variable + equal + exp)
        
        return equation
    
    def _pre(self):
        """
        Do some simple pre-processing on the line
        """
        
        self.line = self.line.replace('>', '')
        self.line = self.line.replace(':=', '=')
        self.line = self.line.replace(';', '')
        self.line = self.line.replace('\n', '')
        
        self.line = self.line.replace('Vector[row]', 'VectorRow')
        
    def parse(self, fname, debug=False, output='parsed.py'):
        
        nr_fails = 0
        result = ''
        with file(fname) as f:
            for i, self.line in enumerate(f.readlines()):
                self._pre()
                try:
                    res_parse = self.bnf.parseString(self.line)
                except ParseException:
                    if debug:
                        print 'xxxx   %s' % (self.line)
                    nr_fails += 1
                res_flat = flatten_to_string(res_parse)
                result += res_flat + '\n'
                if debug:
                    print '%4i   %s' % (i, res_flat)
        
        if debug:
            print
            print '*'*75
            print '%4i   lines where not parsed sucesfully' % nr_fails
        
        with file(output, 'w') as f:
            f.write(result)
            f.flush()

class Tests(unittest.TestCase):
    
    def cases(self):
        ex1 = 'p_omega[1[2]][9a9]'
        ex2 = 'p_omega[1]:=omega[1][1]'
        ex3 = 'R[beta[2]] := Matrix([[cos(beta[2]),0,-sin(beta[2])],[0,1,0],'
        ex3 += '[sin(beta[2]),0,cos(beta[2])]])'
        
        ex4 = '[[cos(beta[2]),0,-sin(beta[2])],[0,1,0],'
        ex4 += '[sin(beta[2]),0,cos(beta[2])]]'
        ex5 = '[cos(beta[2]),0,-sin(beta[2])]'
        
        
        ex61 = 'Stuff(beta[2],0,-beta[2])'
        ex62 = 'Stuff(cos(beta[2]),0,-sin(beta[2]))'
        ex7 = 'p_omega[1],omega[1][1],beta[2]'
        ex8 = 'cos(beta[2]),0,-sin(beta[2])'
        
        ex = 'JT = QQ*m_t*q1*q[1][2][3]'
        res = equation.parseString(ex)
        print res
        print
        print self._flatten_to_string(res)

        ex = 'JT = (1/12)*m_t*(l/2)*(l/(2*qwerty))+9'
        res = equation.parseString(ex)
        print res
        print
        print self._flatten_to_string(res)
        
        # --------------------------------
        ex41 = '[a,b,-3]'
        res41 = matrix.parseString(ex41)
        for k in res41:
            print k
        # --------------------------------
        ex42 = '[ [a,b,-3], [qq,bb,-120] ]'
        res42 = matrix.parseString(ex42)
        for k in res42:
            print k
        # --------------------------------
        
        # --------------------------------        
        arguments_only = OneOrMore( Group(variable + comma) )
        res7 = arguments_only.parseString(ex7)
        for k in res7:
            print k
        # --------------------------------
        ex3 = 'R[beta[2]] = Matrix([[cos(beta[2]),0,-sin(beta[2])],[0,1,0],'
        ex3 += '[sin(beta[2]),0,cos(beta[2])]])'
        res3 = variable.parseString(ex3)
        print res3
        
        # --------------------------------
        ex3 = 'Matrix([[cos(beta[2]),0,-sin(beta[2])],[0,1,0],'
        ex3 += '[sin(beta[2]),0,cos(beta[2])]])'
        res3 = exp.parseString(ex3)
        pprint(res3.asList(), indent=2, depth=20, width=20)
        
        res = exp.parseString(ex62)
        pprint(res.asList(), indent=2, depth=10, width=10)
        # --------------------------------
        
        # --------------------------------
        res61 = command.parseString(ex61)
        for k in res61:
            print k
        # --------------------------------        
        res62 = command.parseString(ex62)
        for k in res62:
            print k
        # --------------------------------
        
        res = equation.parseString(ex2)
        pprint(res.asList(), indent=2, depth=10, width=10)
        
        res = equation.parseString(ex3)
        pprint(res.asList(), indent=2, depth=10, width=10)
        for k in res:
            print k
        
        
        # everything combined
        ex3 = 'dummy = Matrix([[cos(beta[2]),0,-sin(beta[2])],[0,1,0],'
        ex3 += '[sin(beta[2]),VectorAdd(p_omega[1[2]][9a9]*(9+10),beta),cos(beta[2])]])'

if __name__ == '__main__':
    
    fname = 'example2-txtouput.txt'
    mw = MapleGrammer()
    res = mw.parse(fname, debug=True)
    
    

