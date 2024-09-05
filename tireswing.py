import numpy as np
import sympy as sm

A, B, C, k1, k2, x = sm.symbols('A B C k_1 k_2 x', real=True)
M, E, Ix, P, W = sm.symbols('M E I_x P W', real=True)
C1, C2 = sm.symbols('C1 C2')
y = sm.Function('y')(x)
yp = y.diff(x)
ypp = yp.diff(x)
ode = ypp + A*y + B*x + C
# rhs1 = sm.dsolve(ode, y).rhs
# rhs2 = rhs1.subs(sm.sqrt(-A), sm.sqrt(-1)*sm.sqrt(A))
# rhs3 = rhs2.subs(sm.exp(sm.sqrt(-1)*sm.sqrt(A)*x), sm.cos(sm.sqrt(A)*x)+sm.sqrt(-1)*sm.sin(sm.sqrt(A)*x))
# rhs2 = C1*sm.exp(-x*sm.sqrt(-1)*sm.sqrt(A))+C2*sm.exp(x*sm.sqrt(-1)*sm.sqrt(A))-B*x/A - C/A
# rhs3 = C1*(sm.cos(sm.sqrt(A)*x) - sm.sqrt(-1)*sm.sin(sm.sqrt(A)*x)) + C2*(sm.cos(sm.sqrt(A)*x) + sm.sqrt(-1)*sm.sin(sm.sqrt(A)*x)) - B*x/A - C/A

# eqns1 = [rhs1.subs(x,0), rhs1.diff(x).subs(x,0)]
# sols1 = sm.solve(eqns1, C1, C2, dict=True)
# rhs41 = rhs1.subs({C1:sols1[0][C1], C2:sols1[0][C2]})

# eqns3 = [rhs3.subs(x,0), rhs3.diff(x).subs(x,0)]
# sols3 = sm.solve(eqns3, C1, C2, dict=True)
# rhs43 = rhs3.subs({C1:sols3[0][C1], C2:sols3[0][C2]})

# This line combines the three [substitution + solve + back-substitution] steps into one by specifying the initial conditions right at the outset
rhs = sm.dsolve(ode, y, ics={y.subs(x,0):0, yp.subs(x,0):0}).rhs

# Nested set of methods to do the following:
#   Turn sqrt(-A)'s into I*sqrt(A)
#   Turn exp(I*x) into cos(x) + I*sin(x)
#   Simplify
rhs = sm.powdenest(rhs, force=True).rewrite(sm.cos).simplify()



def my_eqns(x, A, B, C, k1, k2):
    y = B/A + k2*sm.sin(sm.sqrt(A)*x) + k1*sm.cos(sm.sqrt(A)*x) + x*C/A
    yp = sm.diff(y, x)
    # yp = -sm.sqrt(A)*k1*sm.sin(sm.sqrt(A)*x) + sm.sqrt(A)*k2*sm.cos(sm.sqrt(A)*x) + C/A
    return y, yp

foo = my_eqns(x, A, B, C, k1, k2)
