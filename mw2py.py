# -*- coding: utf-8 -*-
"""

"""

__author__ = 'David Verelst'
__license__ = 'GPLv3'

from IPython.nbformat import current
import io

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

def flatten_to_string(iterable):
    """
    Flatten the multi-depth iterable and return as a single string
    """
    parsed_line = ''
    for k in flatten(iterable):
        parsed_line += k
    
    return parsed_line

class MapleGrammer:
    
    def __init__(self):
        
        self.bnf = self.BNF()
        self.var_stack = {}
        self.command_stack = {}
        
        self.file_header = """
import numpy as np
import scipy as sp
import sympy as sy
from sympy.matrices import Matrix
from sympy import *
# LaTeX rendered SymPy output
from sympy import init_printing 
"""
        self.maple_2_sympy_mappings = """
def Transpose(matrix):
    return sy.transpose()

def VectorRow(matrix):
    return sy.transpose(sy.Matrix(matrix))

def VectorAdd(a, b):
    return a+b

def Multiply(a, b):
    return a*b
"""
    
    def _add_var(self, strg, loc, toks ):
        """
        Maintain a list of variables so we can define them as proper symbolic
        variables later with SymPy
        """
#        var_stack.append( [f for f in flatten(toks)] )
#        var_stack.append( flatten_to_string(toks) )
        
        varname = flatten_to_string(toks)
        
        # do not include numbers
        try:
            float(varname)
        except ValueError:
            # strip - or + that might be leading a term
            if varname[0] in ['+', '-']:
                self.var_stack[varname[1:]] = strg
            else:
                self.var_stack[varname] = strg
    
    def _add_command(self, strg, loc, toks):
        """
        check which commands have been used in the Maple context
        """
        command = flatten_to_string(toks)
        command = command.split('(')[0]
        if command[0] in ['+', '-']:
            self.command_stack[command[1:]] = strg
        else:
            self.command_stack[command] = strg
        
    
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
        # squary brackets in the variable names are filtered out, so use the
        # l|rsbr_ with the suppress() function
        subname = ZeroOrMore(Group(lsbr_ + variable + rsbr_))
        variable << Group( (varname|nums) + subname )
        # keep track of all the variables added
        variable.setParseAction( self._add_var )
        
        # same with commands, matrices: they can be nested.
        # mind the order: first command, than variable. The variable definition
        # will be satisfied just before reaching the opening bracket and that
        # will result in a false positive for the variable
        args = OneOrMore( Group((matrix|command|mixed_terms|variable)+comma) )
        command << Group( (varname + Group(lpar + args + rpar)) )
        command.setParseAction(self._add_command)
        
        # matrices are also nested structures
        matrix << Group(lsbr + args + rsbr)
        
        # terms in parenthesis can hold other nested structures
        terms = OneOrMore( Group( (variable|terms_par) + op_opt) )
        terms_par << OneOrMore( Group(lpar + terms + rpar + op_opt) )
        mixed_terms << OneOrMore( (terms_par|terms) + op_opt )
        
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
    
    def define_sympy_symbols(self):
        
        # define all the sympy variables
        self.sympy_vars = ''
        for var in sorted(self.var_stack.keys()):
            self.sympy_vars += "%s = Symbol('%s')\n" % (var, var)
#        result.append(["%s = Symbol('%s')\n" % var for var in self.var_stack])
    
    def parse_line(self, line):
        """
        """
        
        self.line = line
        self._pre()
        return self.bnf.parseString(self.line)
    
    def parse_file(self, fname, debug=False):
        
        nr_fails = 0
        result = [self.file_header]
        result.append(self.maple_2_sympy_mappings)
        # empty cell to add sympy vars later
        result.append('')
        with file(fname) as f:
            for i, self.line in enumerate(f.readlines()):
                try:
                    res_parse = self.parse_line(self.line)
                except ParseException:
                    res_parse = '#' + self.line
                    nr_fails += 1
                res_flat = flatten_to_string(res_parse)
                result.append(res_flat)
                if debug:
                    print '%4i   %s' % (i, res_flat)
        
        self.define_sympy_symbols()
        result[2] = self.sympy_vars
        
        if debug:
            print
            print '*'*75
            print '%4i   lines where not parsed sucesfully' % nr_fails

        return result
        
    def write_py(self, result, output='parsed-output.py'):
        
        with file(output, 'w') as f:
            # since writelines doesn't add newlines, add them manually
            f.writelines(map(lambda line: line+'\n', result))
            f.flush()
    
    def write_notebook(self, result, output='parsed-output.ipynb'):
        """
        based on StackOverflow quetsion and answer:
            http://stackoverflow.com/questions/17972273/
            is-there-a-ipython-notebook-api
        """
        
        # open an empty notebook
        notebook = current.reads('', format='py')
        # add all elements as seperate cells
        notebook['worksheets'][0]['cells'] = list(map(current.new_code_cell,
                                                      result))
        # Save the notebook
        with io.open(output, 'w', encoding='utf-8') as f:
            current.write(notebook, f, format='ipynb')
    

class Tests(unittest.TestCase):
    
    def setUp(self):
        self.mw = MapleGrammer()
        
    
    def test_case1(self):
        line = '> omega[beta[2]] := Vector[row]([0,beta[flux[2]],0]);'
        result = flatten_to_string(self.mw.parse_line(line))
        truth = 'omegabeta2=VectorRow([0,betaflux2,0])'
        self.assertEqual(result, truth)
        
        line = '> r_0 := Vector[row]([l/2,0,0]);'
        result = flatten_to_string(self.mw.parse_line(line))
        truth = 'r_0=VectorRow([l/2,0,0])'
        self.assertEqual(result, truth)
        
#> r_cg_beta := Vector[row]([R/2,0,0]);
#> r_cg_T := Vector[row]([l/4,0,0]);
#> r_3_beta_1 := Vector[row]([(-e*cos(psi)),(e*sin(psi)),0]);
#> r_3_beta_2 := Vector[row]([(e*cos(psi)),(-e*sin(psi)),0]);
#> r_3_beta_1_flux := Vector[row]([(e*sin(psi)*psi[flux]),(e*cos(psi)*psi[flux]),0]);
    
    def test_case2(self):
        line = '> dT_dbeta1_t := subs({beta[flux[1]]=beta[flux[1]](t),'
        line += 'beta[flux[2]]=beta[flux[2]](t),rho[flux]=rho[flux](t), '
        line += 'psi[flux]=Omega, beta[1]=beta[1](t), beta[2]=beta[2](t), '
        line += 'rho=rho(t), psi=Omega*t},dT_dbeta1);'
    
    def test_case3(self):
        line = '> omega[rho[1]] := Multiply(omega[rho],Multiply('
        line += 'Transpose(R[psi[1]]),Transpose(R[beta[1]])));'
        
        truth = 'omegarho1=Multiply(omegarho,Multiply(Transpose(Rpsi1),'
        truth += 'Transpose(Rbeta1)))'
    
    def test_case4(self):
        """
        """
        
        line = '> R_cg_blade_1_0 := r_0 + Multiply(r_2,R[rho])'
        line += '+ Multiply(r_3_beta_1,R[rho])'
        line += '+ Multiply(r_cg_beta,Multiply(R[beta[1]],'
        line += 'Multiply(R[psi[1]],R[rho])));'
        
        truth = 'R_cg_blade_1_0=r_0+Multiply'
    
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
    
#    unittest.main()
    
    mw = MapleGrammer()
    
    # single line parsing
    line = '> r_0 := Vector[row]([l/2,0,0]);'
    line = '> omega[rho[1]] := Multiply(omega[rho],Multiply('
    line += 'Transpose(R[psi[1]]),Transpose(R[beta[1]])));'
    result = mw.parse_line(line)
    res_flat = flatten_to_string(result)
    print line
    print res_flat
    
    # parse a complete file
    fname = 'example2-txtouput.txt'
    res = mw.parse_file(fname, debug=True)
    mw.write_py(res)
    mw.write_notebook(res)
    
    print 
    print 'variables used'
    print '-'*22
    for var in sorted(mw.var_stack.keys()): print var
    
    print
    print 'commands used by Maple'
    print '-'*22
    for k in sorted(mw.command_stack.keys()): print k
    
    

