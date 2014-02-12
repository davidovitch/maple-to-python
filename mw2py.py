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
        self.vector_stack = {}
        self.command_stack = {}
        self.var_replace = []
        self.correct_var_index = False
        
        self.file_header = """
import numpy as np
import scipy as sp
import sympy as sy
from sympy.matrices import Matrix
from sympy import *
# LaTeX rendered SymPy output
from sympy import init_printing
init_printing()
"""
        self.maple_2_sympy_mappings = """
def Transpose(matrix):
    return sy.transpose(matrix)

def VectorRow(vector_list):
    return sy.transpose(Matrix(vector_list))

def VectorAdd(a, b):
    return a+b

def Multiply(a, b):
    return a*b
"""
   
    def _add_vector(self, strg, loc, toks):
        """
        """
        varname = flatten_to_string(toks)
        self.vector_stack[varname] = strg
   
    def _add_var(self, strg, loc, toks):
        """
        Maintain a list of variables so we can define them as proper symbolic
        variables later with SymPy
        """
#        var_stack.append( [f for f in flatten(toks)] )
#        var_stack.append( flatten_to_string(toks) )
        
        varname = flatten_to_string(toks)
        
        # if the last character is and integer, and if the base name has been
        # declared before, we have have to replace the var name
        try:
            i = int(varname[-1])
            if varname[:-1] in self.var_stack:
                # do not add the variable, and mark for replacement
                newname = '%s[%i]' % (varname[:-1], i-1)
                print '+++++++++', varname, '---->', newname
                self.var_replace.append( [varname, newname])
                self.correct_var_index = True
                return
        except ValueError:
            self.correct_var_index = False
        
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
        # Should be optional for 1x [] group. With ZeroOrMore it will also 
        # match omega[1][1], which is actually pointing to an index
        # FIXME: the solution should be: if the variable xyz[subscript] has
        # been defined earlier in the script, then xyz[subscript][i] is and
        # index and NOT part of the variable name. This means that we should
        # build up the entire matching process through a stack variable
        # that is callend upon setParseAction.
        # Note that also omegarho3 should be omegarho[2] and D_tower2 as well
#        subname = ZeroOrMore( Group(lsbr_ + variable + rsbr_) )
#        variable << Group( (varname|nums) + subname )
        # subname optional because with ZeroOrMore you get the index:
        # qwerty[sub[subsub]][index]
        subname = Group(lsbr_ + variable + rsbr_)
        variable << Group( (varname|nums) + Optional(subname) )
        # keep track of all the variables added
        variable.setParseAction( self._add_var )
#        self.variable = variable
        
        # a vector is actually a special case of a variable. The problem is
        # that most vectors are also valid variables. However, one case should
        # only result in a vector: qwerty[1][1]
        vector = Group(variable + lsbr + nums + rsbr).setParseAction( self._add_vector )
        
        # same with commands, matrices: they can be nested.
        # mind the order: first command, than variable. The variable definition
        # will be satisfied just before reaching the opening bracket and that
        # will result in a false positive for the variable
        args = OneOrMore( Group((matrix|command|mixed_terms|variable|vector)+comma) )
        command << Group( (varname + Group(lpar + args + rpar)) )
        command.setParseAction(self._add_command)
        
        # matrices are also nested structures
        matrix << Group(lsbr + args + rsbr)
        
        # terms in parenthesis can hold other nested structures
        # FIXME: variable here shouldn't be required?
        terms = OneOrMore( Group( (command|variable|vector|terms_par) + op_opt) )
        terms_par << OneOrMore( Group(lpar + terms + rpar + op_opt) )
        mixed_terms << OneOrMore( (terms_par|terms) + op_opt )
        
        # an expression of the form
        exp = (command | mixed_terms)
        # and finally we can assemble an equation
        equation = ( (variable|vector) + equal + exp)
        
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
        
#        # parse the LHS first, that can only be a variable
#        try:
#            lhs, rhs = self.line.split('=')
#            lhs_parse = self.variable.parseString(lhs)
#            lhs_parse = flatten_to_string(lhs_parse)
#            
#            rhs_parse = self.variable.parseString(lhs)
#            rhs_parse = flatten_to_string(rhs_parse)
#            
#        except ValueError:
        
        return self.bnf.parseString(self.line)
    
    def monkey_patch(self, line):
        """
        check if the basename has beend delcared before
        """
        
        if self.correct_var_index:
            target = min(len(self.var_replace), 12)
            for i in range(-1, -target, -1):
                line = line.replace(self.var_replace[i][0], 
                                    self.var_replace[i][1])
            return line
        else:
            return line
    
    def parse_file(self, fname, debug=False):
        
        nr_fails = 0
        result = [self.file_header]
        result.append(self.maple_2_sympy_mappings)
        # empty cell to add sympy vars later
        result.append('')
        with file(fname) as f:
            for i, self.line in enumerate(f.readlines()):
                # lines to ignore, dirty hacks
                if len(self.line) > 1000:
                    res_parse = '# too long...'
                elif self.line.startswith('Equation'):
                    continue
                elif self.line.find(':=') == -1:
                    res_parse = '#' + self.line
                else:
                    try:
                        res_parse = self.parse_line(self.line)
                    except ParseException:
                        res_parse = '#' + self.line
                        nr_fails += 1
                res_flat = flatten_to_string(res_parse)
                res_flat = self.monkey_patch(res_flat)
                result.append(res_flat)
                if debug:
                    print '%4i   %s' % (i, res_flat)
        
        # replace the variables that have their index included
        
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
    
    def print_test_info(self):
        """
        """
        print 
        print 'variables used'
        print '-'*22
        for var in sorted(self.mw.var_stack.keys()): 
            print var
        
        print 
        print 'vectors found'
        print '-'*22
        for var in sorted(self.mw.vector_stack.keys()): 
            print var
        
        print
        print 'commands used by Maple'
        print '-'*22
        for k in sorted(self.mw.command_stack.keys()): 
            print k
    
    def test_case1(self):
        """
        variables with multi depth subscripts
        """
        line = '> omega[beta[2]] := Vector[row]([0,beta[flux[2]],0]);'
        truth = 'omegabeta2=VectorRow([0,betaflux2,0])'
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
        
        line = '> r_0 := Vector[row]([l/2,0,0]);'
        truth = 'r_0=VectorRow([l/2,0,0])'
        result = flatten_to_string(self.mw.parse_line(line))
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
        """
        Multi depth commands
        """
        line = '> omega[rho[1]] := Multiply(omega[rho],Multiply('
        line += 'Transpose(R[psi[1]]),Transpose(R[beta[1]])));'
        
        truth = 'omegarho1=Multiply(omegarho,Multiply(Transpose(Rpsi1),'
        truth += 'Transpose(Rbeta1)))'
        
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
    
    def test_case4(self):
        """
        """
        
        line = '> R_cg_blade_1_0 := r_0 + Multiply(r_2,R[rho])'
        line += '+ Multiply(r_3_beta_1,R[rho])'
        line += '+ Multiply(r_cg_beta,Multiply(R[beta[1]],'
        line += 'Multiply(R[psi[1]],R[rho])));'
        
        truth = 'R_cg_blade_1_0=r_0+Multiply(r_2,Rrho)'
        truth += '+Multiply(r_3_beta_1,Rrho)'
        truth += '+Multiply(r_cg_beta,Multiply(Rbeta1,'
        truth += 'Multiply(Rpsi1,Rrho)))'
        
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
    
    def test_case5(self):
        line = '> JT := (1/12)*m_t*(l/2)*(l/2);'
        truth = 'JT=(1/12)*m_t*(l/2)*(l/2)'
        
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
        
        line = 'JT := (1/12)*m_t*(l/2)*(l/(2*qwerty))+9;'
        truth = 'JT=(1/12)*m_t*(l/2)*(l/(2*qwerty))+9'
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
    
    def test_case6(self):
        """
        """
        line = '> omega[beta[1]] := Vector[row]([0,beta[flux[1]],0]);'
        truth = 'omegabeta1=VectorRow([0,betaflux1,0])'
        
        result = flatten_to_string(self.mw.parse_line(line))
        self.assertEqual(result, truth)
        self.assertTrue('betaflux1' in self.mw.var_stack.keys())
        variables = self.mw.var_stack.keys()
        try:
            self.assertFalse('beta1' in variables)
            self.assertFalse('flux1' in variables)
        except AssertionError as e:
            self.print_test_info()
            raise e
    
    def test_case7(self):
        pass
    
    def test_case8(self):
        """
        find vector_index, multiple recursive terms
        """
        
        # very simple case
        line = '> T_rot_blade_2[a][1] := omega[2][1]'#+'
        truth  = 'T_rot_blade_2a[1]=omega2[1]'#+'
        result = flatten_to_string(self.mw.parse_line(line))
        self.print_test_info()
        self.assertEqual(result, truth)
        
        # medium case
        line = '> T_rot_blade_2[a][1] := omega[2][1]+'
        line += '(D_blade_2[2]*omega[2][2])+(D_blade_2[3]*omega[2][3]);'
        truth  = 'T_rot_blade_2a[1]=omega2[1]+'
        truth += '(D_blade_2[2]*omega2[2])+(D_blade_2[3]*omega[2][3]);'
        result = flatten_to_string(self.mw.parse_line(line))
        self.print_test_info()
        self.assertEqual(result, truth)
        
        line = '> T_rot_blade_2 := 0.5*((D_blade_2[1]*omega[2][1])+'
        line += '(D_blade_2[2]*omega[2][2])+(D_blade_2[3]*omega[2][3]));'
        truth  = 'T_rot_blade_2=0.5*((D_blade_2[1]*omega2[1])+'
        truth += '(D_blade_2[2]*omega2[2])+(D_blade_2[3]*omega[2][3]));'
        result = flatten_to_string(self.mw.parse_line(line))
        self.print_test_info()
        self.assertEqual(result, truth)
        
    
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
    
    unittest.main()
    
    mw = MapleGrammer()
    
    # single line parsing
    line = '> r_0 := Vector[row]([l/2,0,0]);'
    line = '> omega[rho[1]] := Multiply(omega[rho],Multiply('
    line += 'Transpose(R[psi[1]]),Transpose(R[beta[1]])));'
    line = '> omega[beta[1]] := Vector[row]([0,beta[flux[1]],0]);'
    result = mw.parse_line(line)
    res_flat = flatten_to_string(result)
    print line
    print res_flat
    
    print 
    print 'variables used'
    print '-'*22
    for var in sorted(mw.var_stack.keys()): print var
    
    print
    print 'commands used by Maple'
    print '-'*22
    for k in sorted(mw.command_stack.keys()): print k
    
    # parse a complete file
    fname = 'example2-txtouput-full.txt'
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
    
    unittest.main()

