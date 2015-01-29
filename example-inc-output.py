
import numpy as np
import scipy as sp
import sympy as sy
from sympy.matrices import Matrix
from sympy import *
# LaTeX rendered SymPy output
from sympy import init_printing
init_printing()


def Transpose(matrix):
    return sy.transpose(matrix)

def VectorRow(vector_list):
    return sy.transpose(Matrix(vector_list))

def VectorAdd(a, b):
    return a+b

def Multiply(a, b):
    return a*b

D_blade_1 = Symbol('D_blade_1')
D_blade_2 = Symbol('D_blade_2')
D_tower = Symbol('D_tower')
JB = Symbol('JB')
JT = Symbol('JT')
Jblade = Symbol('Jblade')
Jtower = Symbol('Jtower')
Omega = Symbol('Omega')
Omega1 = Symbol('Omega1')
OmegaXomega1 = Symbol('OmegaXomega1')
OmegaXomega2 = Symbol('OmegaXomega2')
OmegaXomegarho = Symbol('OmegaXomegarho')
R = Symbol('R')
R_cg_T = Symbol('R_cg_T')
R_cg_T_0 = Symbol('R_cg_T_0')
R_cg_T_flux = Symbol('R_cg_T_flux')
R_cg_beta_1 = Symbol('R_cg_beta_1')
R_cg_beta_1_flux = Symbol('R_cg_beta_1_flux')
R_cg_beta_2 = Symbol('R_cg_beta_2')
R_cg_beta_2_flux = Symbol('R_cg_beta_2_flux')
R_cg_blade_1_0 = Symbol('R_cg_blade_1_0')
R_cg_blade_2_0 = Symbol('R_cg_blade_2_0')
Rbeta1 = Symbol('Rbeta1')
Rbeta2 = Symbol('Rbeta2')
Rpsi1 = Symbol('Rpsi1')
Rpsi2 = Symbol('Rpsi2')
Rrho = Symbol('Rrho')
T1_eq_beta1full = Symbol('T1_eq_beta1full')
T1_eq_beta2full = Symbol('T1_eq_beta2full')
T1_eq_rhofull = Symbol('T1_eq_rhofull')
T2_eq_beta1full = Symbol('T2_eq_beta1full')
T2_eq_beta2full = Symbol('T2_eq_beta2full')
T2_eq_rhofull = Symbol('T2_eq_rhofull')
T3_eq_beta1full = Symbol('T3_eq_beta1full')
T3_eq_beta2full = Symbol('T3_eq_beta2full')
T3_eq_rhofull = Symbol('T3_eq_rhofull')
T4_eq_rhofull = Symbol('T4_eq_rhofull')
T5_eq_rhofull = Symbol('T5_eq_rhofull')
T_blade_1 = Symbol('T_blade_1')
T_blade_2 = Symbol('T_blade_2')
T_rot_blade_1 = Symbol('T_rot_blade_1')
T_rot_blade_2 = Symbol('T_rot_blade_2')
T_rot_tower = Symbol('T_rot_tower')
T_total = Symbol('T_total')
T_tower = Symbol('T_tower')
V_g = Symbol('V_g')
V_springs = Symbol('V_springs')
V_total = Symbol('V_total')
Xomega1 = Symbol('Xomega1')
Xomega2 = Symbol('Xomega2')
Xomegarho = Symbol('Xomegarho')
aa_FF = Symbol('aa_FF')
bb_FF_1 = Symbol('bb_FF_1')
bb_FF_2 = Symbol('bb_FF_2')
beta1 = Symbol('beta1')
beta2 = Symbol('beta2')
betaflux1 = Symbol('betaflux1')
betaflux2 = Symbol('betaflux2')
blade = Symbol('blade')
dT_dbeta1 = Symbol('dT_dbeta1')
dT_dbeta1_t = Symbol('dT_dbeta1_t')
dT_dbeta1_tt = Symbol('dT_dbeta1_tt')
dT_dbeta1flux = Symbol('dT_dbeta1flux')
dT_dbeta1flux_t = Symbol('dT_dbeta1flux_t')
dT_dbeta1flux_tt = Symbol('dT_dbeta1flux_tt')
dT_dbeta2 = Symbol('dT_dbeta2')
dT_dbeta2_t = Symbol('dT_dbeta2_t')
dT_dbeta2_tt = Symbol('dT_dbeta2_tt')
dT_dbeta2flux = Symbol('dT_dbeta2flux')
dT_dbeta2flux_t = Symbol('dT_dbeta2flux_t')
dT_dbeta2flux_tt = Symbol('dT_dbeta2flux_tt')
dT_drho = Symbol('dT_drho')
dT_drho_t = Symbol('dT_drho_t')
dT_drho_tt = Symbol('dT_drho_tt')
dT_drhoflux = Symbol('dT_drhoflux')
dT_drhoflux_t = Symbol('dT_drhoflux_t')
dT_drhoflux_tt = Symbol('dT_drhoflux_tt')
dV_dbeta1 = Symbol('dV_dbeta1')
dV_dbeta1_t = Symbol('dV_dbeta1_t')
dV_dbeta1_tt = Symbol('dV_dbeta1_tt')
dV_dbeta2 = Symbol('dV_dbeta2')
dV_dbeta2_t = Symbol('dV_dbeta2_t')
dV_dbeta2_tt = Symbol('dV_dbeta2_tt')
dV_drho = Symbol('dV_drho')
dV_drho_t = Symbol('dV_drho_t')
dV_drho_tt = Symbol('dV_drho_tt')
ddt_dT_dbeta1flux = Symbol('ddt_dT_dbeta1flux')
ddt_dT_dbeta2flux = Symbol('ddt_dT_dbeta2flux')
ddt_dT_drhoflux = Symbol('ddt_dT_drhoflux')
diff = Symbol('diff')
e = Symbol('e')
eq_beta1full = Symbol('eq_beta1full')
eq_beta2full = Symbol('eq_beta2full')
eq_mo_beta1_FF_lin = Symbol('eq_mo_beta1_FF_lin')
eq_mo_beta2_FF_lin = Symbol('eq_mo_beta2_FF_lin')
eq_mo_rho_FF_lin = Symbol('eq_mo_rho_FF_lin')
eq_rhofull = Symbol('eq_rhofull')
flux = Symbol('flux')
flux1 = Symbol('flux1')
full = Symbol('full')
g = Symbol('g')
g_ = Symbol('g_')
horner = Symbol('horner')
ka = Symbol('ka')
kb = Symbol('kb')
l = Symbol('l')
m_b = Symbol('m_b')
m_t = Symbol('m_t')
omega1 = Symbol('omega1')
omega2 = Symbol('omega2')
omegaOmega1 = Symbol('omegaOmega1')
omegaOmega2 = Symbol('omegaOmega2')
omegabeta1 = Symbol('omegabeta1')
omegabeta2 = Symbol('omegabeta2')
omegarho = Symbol('omegarho')
p_omega1 = Symbol('p_omega1')
p_omega2 = Symbol('p_omega2')
psi = Symbol('psi')
psi1 = Symbol('psi1')
psiflux = Symbol('psiflux')
q_omega1 = Symbol('q_omega1')
q_omega2 = Symbol('q_omega2')
r_0 = Symbol('r_0')
r_2 = Symbol('r_2')
r_3_beta_1 = Symbol('r_3_beta_1')
r_3_beta_1_flux = Symbol('r_3_beta_1_flux')
r_3_beta_2 = Symbol('r_3_beta_2')
r_3_beta_2_flux = Symbol('r_3_beta_2_flux')
r_cg_T = Symbol('r_cg_T')
r_cg_beta = Symbol('r_cg_beta')
r_omega1 = Symbol('r_omega1')
r_omega2 = Symbol('r_omega2')
rho = Symbol('rho')
rhoflux = Symbol('rhoflux')
solve = Symbol('solve')
subs = Symbol('subs')
t = Symbol('t')
tower = Symbol('tower')

#load linear algebra package

#> restart;

#> with(CodeGeneration):

#> with(LinearAlgebra):

#Transformation matrices for tower middle

Rrho=Matrix([[cos(rho),0,-sin(rho)],[0,1,0],[sin(rho),0,cos(rho)]])
#Transformation matrices for blade 1

Rpsi1=Matrix([[-cos(psi),sin(psi),0],[-sin(psi),-cos(psi),0],[0,0,1]])
Rbeta1=Matrix([[cos(beta1),0,-sin(beta1)],[0,1,0],[sin(beta1),0,cos(beta1)]])
#Transformation matrices for blade 2

Rpsi2=Matrix([[cos(psi),-sin(psi),0],[sin(psi),cos(psi),0],[0,0,1]])
Rbeta2=Matrix([[cos(beta2),0,-sin(beta2)],[0,1,0],[sin(beta2),0,cos(beta2)]])
#omega of blade 1, in function of {E_beta_1}

#omega of the tower, in function of {E_1}

omegarho=VectorRow([0,rhoflux,0])
omegarho[0]=Multiply(omegarho,Multiply(Transpose(Rpsi[0]),Transpose(Rbeta1)))
omegaOmega1=Multiply(VectorRow([0,0,-Omega]),Transpose(Rbeta1))
omegabeta1=VectorRow([0,betaflux[0],0])
omega1=VectorAdd(VectorAdd(omegarho[0],omegaOmega[0]),omegabeta1)
p_omega1=omega1
q_omega1=omega1
r_omega1=omega1
#omega of blade 2, in function of {E_beta_2}

omegarho[1]=Multiply(omegarho,Multiply(Transpose(Rpsi[1]),Transpose(Rbeta2)))
omegaOmega2=Multiply(VectorRow([0,0,-Omega]),Transpose(Rbeta2))
omegabeta2=VectorRow([0,betaflux[1],0])
omega2=VectorAdd(VectorAdd(omegarho[1],omegaOmega[1]),omegabeta2)
p_omega2=omega2
q_omega2=omega2
r_omega2=omega2
#rotation operators

OmegaXomegarho=Matrix([[0,0,-rhoflux],[0,0,0],[rhoflux,0,0]])
OmegaXomega1=Matrix([[0,r_omega1,-q_omega1],[-r_omega1,0,p_omega1],[q_omega1,-p_omega1,0]])
OmegaXomega2=Matrix([[0,r_omega2,-q_omega2],[-r_omega2,0,p_omega2],[q_omega2,-p_omega2,0]])
#Position vectors for the cg of the blades and the cg of the 2nd tower section

#All captial R position vectors here are expressed in {E_1}

#all r position vectors are in local reference frame

r_2=VectorRow([l/2,0,0])
r_0=VectorRow([l/2,0,0])
r_cg_beta=VectorRow([R/2,0,0])
r_cg_T=VectorRow([l/4,0,0])
r_3_beta_1=VectorRow([(-e*cos(psi)),(e*sin(psi)),0])
r_3_beta_2=VectorRow([(e*cos(psi)),(-e*sin(psi)),0])
r_3_beta_1_flux=VectorRow([(e*sin(psi)*psiflux),(e*cos(psi)*psiflux),0])
r_3_beta_2_flux=VectorRow([(-e*sin(psi)*psiflux),(-e*cos(psi)*psiflux),0])
R_cg_beta_1=Multiply(r_0,Transpose(Rrho))
R_cg_beta_2=Multiply(r_0,Transpose(Rrho))
R_cg_T=Multiply(r_0,Transpose(Rrho))
R_cg_beta_1_flux=Multiply(r_2,OmegaXomegarho)
R_cg_beta_2_flux=Multiply(r_2,OmegaXomegarho)
R_cg_T_flux=Multiply(r_cg_T,OmegaXomegarho)
#Inertia matrices of cg

#> JB = (1/12)*m_b*R*R;

Jblade=Matrix([[0,0,0],[0,JB,0],[0,0,JB]])
#> JT = (1/12)*m_t*(l/2)*(l/2);

Jtower=Matrix([[0,0,0],[0,JT,0],[0,0,JT]])
#Kinetic energy: translational part

T_blade_1=0.5*m_b*((R_cg_beta_1_flux[0]*R_cg_beta_1_flux[0])+(R_cg_beta_1_flux[1]*R_cg_beta_1_flux[1])+(R_cg_beta_1_flux[2]*R_cg_beta_1_flux[2]))
T_blade_2=0.5*m_b*((R_cg_beta_2_flux[0]*R_cg_beta_2_flux[0])+(R_cg_beta_2_flux[1]*R_cg_beta_2_flux[1])+(R_cg_beta_2_flux[2]*R_cg_beta_2_flux[2]))
T_tower=0.5*m_t*((R_cg_T_flux[0]*R_cg_T_flux[0])+(R_cg_T_flux[1]*R_cg_T_flux[1])+(R_cg_T_flux[2]*R_cg_T_flux[2]))
#Kinetic energy: rotational part

D_blade_1=Multiply(omega1,Jblade)
D_blade_2=Multiply(omega2,Jblade)
D_tower=Multiply(omegarho,Jtower)
T_rot_blade_1=0.5*
T_rot_blade_2=0.5*
T_rot_tower=0.5*
#Total kinetic energy

T_total=T_blade_1+T_blade_2+T_tower+T_rot_blade_1+T_rot_blade_2+T_rot_tower
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

#Potential energy

g=VectorRow([g_,0,0])
#Position vectors expressed in {E_0} for the potential energy

R_cg_blade_1_0=r_0+Multiply(r_2,Rrho)+Multiply(r_3_beta_1,Rrho)+Multiply(r_cg_beta,Multiply(Rbeta1,Multiply(Rpsi1,Rrho)))
R_cg_blade_2_0=r_0+Multiply(r_2,Rrho)+Multiply(r_3_beta_2,Rrho)+Multiply(r_cg_beta,Multiply(Rbeta2,Multiply(Rpsi2,Rrho)))
R_cg_T_0=r_0+Multiply(r_cg_T,Rrho)
V_g=(m_b*g[0]*(R_cg_blade_1_0[0]+R_cg_blade_2_0[0]))+(m_t*g[0]*R_cg_T_0[0])
V_springs=1/2*kb*beta1*beta1+1/2*kb*beta2*beta2+1/2*ka*rho*rho
V_total=V_g+V_springs
#> `starting length` = length(V_total), `converting to horner`=length(convert(V_total,horner)),  factoring=length(factor(V_total)), `simplify length` = length(symplify(V_total));

#Full equations of motion using the Lagrangian

#eq of motion for beta_1

dT_dbeta1flux=diff(T_total,betaflux[0])
dT_dbeta1=diff(T_total,beta1)
dV_dbeta1=diff(V_total,beta1)
#indicate time dependencies

#indicate time dependencies

dT_dbeta1_t=subs
dV_dbeta1_t=subs
dT_dbeta1flux_t=subs
#change q_flux to d/dt(q)

dT_dbeta1_tt=subs
dV_dbeta1_tt=subs
dT_dbeta1flux_tt=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

dT_dbeta1flux_tt=convert(dT_dbeta1flux_tt,horner)
#differentiate with respect to t and construct equations of motion for beta 1

ddt_dT_dbeta1flux=diff(dT_dbeta1flux_tt,t)
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

ddt_dT_dbeta1flux=convert(ddt_dT_dbeta1flux,horner)
eq_beta1full=ddt_dT_dbeta1flux-dT_dbeta1_tt+dV_dbeta1_tt
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

eq_beta1full=convert(eq_beta1full,horner)
eq_beta1full=solve(eq_beta1full,diff(beta1(t),t,t))
#simplify expression and lose time dependencie notation and d/dt

T1_eq_beta1full=subs
T2_eq_beta1full=subs
T3_eq_beta1full=subs
bb_FF_1=T3_eq_beta1full
#eq of motion for beta_2

dT_dbeta2flux=diff(T_total,betaflux[1])
dT_dbeta2=diff(T_total,beta2)
dV_dbeta2=diff(V_total,beta2)
#indicate time dependencies

#indicate time dependencies

dT_dbeta2_t=subs
dV_dbeta2_t=subs
dT_dbeta2flux_t=subs
#change q_flux to d/dt(q)

dT_dbeta2_tt=subs
dV_dbeta2_tt=subs
dT_dbeta2flux_tt=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

dT_dbeta2flux_tt=convert(dT_dbeta2flux_tt,horner)
#differentiate with respect to t and construct equations of motion for beta 2

ddt_dT_dbeta2flux=diff(dT_dbeta2flux_tt,t)
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

ddt_dT_dbeta2flux=convert(ddt_dT_dbeta2flux,horner)
eq_beta2full=ddt_dT_dbeta2flux-dT_dbeta2_tt+dV_dbeta2_tt
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

eq_beta2full=convert(eq_beta2full,horner)
eq_beta2full=solve(eq_beta2full,diff(beta2(t),t,t))
#simplify expression and lose time dependencie notation and d/dt

T1_eq_beta2full=subs
T2_eq_beta2full=subs
T3_eq_beta2full=subs
bb_FF_2=T3_eq_beta2full
#eq of motion for rho

dT_drhoflux=diff(T_total,rhoflux)
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

#replace first rho[flux] with rho_f, otherwise maple will display the partial derivative of d/drho (rho[flux]) which is ofcourse zero, but maple doesn't know

dT_drho=diff
dV_drho=diff(V_total,rho)
#indicate time dependencies

#indicate time dependencies

dT_drho_t=subs
dV_drho_t=subs
dT_drhoflux_t=subs
#change q_flux to d/dt(q)

dT_drho_tt=subs
dV_drho_tt=subs
dT_drhoflux_tt=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

#differentiate with respect to t and construct equations of motion for rho

ddt_dT_drhoflux=diff(dT_drhoflux_tt,t)
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

eq_rhofull=ddt_dT_drhoflux-dT_drho_tt+dV_drho_tt
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

eq_rhofull=solve(eq_rhofull,diff(rho(t),t,t))
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

eq_rhofull=convert(eq_rhofull,horner)
#simplify expression for matlab en lose time dependencie notation and d/dt

T1_eq_rhofull=subs
#> length(%);

T2_eq_rhofull=subs
T3_eq_rhofull=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

T3_eq_rhofull=collect(T3_eq_rhofull,aa_FF)
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

#Solve again for aa_FF: since bb_FF_1 and bb_FF_2 are now substituted and are also functions aa_FF - CAN NOT SIMPLIFY FURTHER

T4_eq_rhofull=solve
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

T5_eq_rhofull=convert(T4_eq_rhofull,horner)
#> length(%);

#Linearized equations of motion using the Lagrangian 

#lin eq of motion for beta_1

eq_mo_beta1_FF_lin=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%)), `simplify length` = length(symplify(%));

#lin eq of motion for beta_2

eq_mo_beta2_FF_lin=subs
#> `starting length` = length(%), `converting to horner`=length(convert(%,horner)),  factoring=length(factor(%));

#lin eq of motion for rho

eq_mo_rho_FF_lin=subs
#Matlab code

#rho

#> Matlab(T5_eq_rho[full], resultname="aa_FF");

# too long...
#> Matlab(eq_mo_rho_FF_lin, resultname="aa_FF_lin");

# too long...
#beta_1

#> Matlab(bb_FF_1, resultname="bb_FF_1");

# too long...
#> Matlab(eq_mo_beta1_FF_lin, resultname="bb_FF_1_lin");

# too long...
#beta_2

#> Matlab(bb_FF_2, resultname="bb_FF_2");

# too long...
#> Matlab(eq_mo_beta2_FF_lin, resultname="bb_FF_2_lin");

# too long...
#

