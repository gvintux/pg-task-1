from sympy import *

from common.models.kelvin_fought.deflection.symbolic import deflection_solve as defsolve

D = Symbol('D', real=True, nonzero=True)
pi = Symbol('rho_i', real=True, nonzero=True)
pw = Symbol('rho_w', real=True, nonzero=True)
h = Symbol('h', real=True, nonzero=True)
g = Symbol('g', real=True, nonzero=True)
t = Symbol('t', real=True, nonzero=True)
q = Symbol('q', real=True, nonzero=True)
r = Symbol('r', real=True, nonzero=True)
mu = Symbol('mu', real=True, nonzero=True)
nu = Symbol('nu', real=True, nonzero=True)
H = Symbol('H', real=True, nonzero=True)
x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
v = Symbol('v', real=True)
tf = Symbol('tau_phi', real=True, nonzero=True)
lm = Symbol('lambda', real=True, nonzero=True)
eta = Symbol('eta', real=True, nonzero=True)
delta = Symbol('delta', real=True)
xi = Symbol('xi', real=True, nonzero=True)
chi = Symbol('chi', real=True, nonzero=True)
a = Symbol('a', real=True, nonzero=True)
b = Symbol('b', real=True, nonzero=True)
A = Symbol('A', real=True, nonzero=True)
B = Symbol('B', real=True, nonzero=True)
dPhidt = Symbol('dPhi_dt', real=True, nonzero=True)
P = Symbol('P', real=True, nonzero=True)
k = Symbol('k')

w1 = Function('w')(x, y)
w2 = Function('w')(t)
w = w1 + w2
w = Function('w')(x, y, t)
Phi1 = Function('Phi')(x, y, z)
Phi2 = Function('Phi')(t)
Phi = Phi1 + Phi2
Phi = Function('Phi')(x, y, z, t)

def tension_solve(**specs):
    if 'lm' in specs:
        global lm
        lm = 0
    if 'eta' in specs:
        global eta
        eta = 0
    w, Aa, Bb = defsolve()
    d2w_dx = diff(w, x, 2)
    d2w_dy = diff(w, y, 2)
    d2w_dx_dt = diff(d2w_dx, t)
    d2w_dy_dt = diff(d2w_dy, t)
    d2w_dx_dy = diff(w, x, 1, y, 1)
    d2w_dx_dy_dt = diff(d2w_dx_dy, t)

    Mx = -D * (d2w_dx + mu * d2w_dy + tf * (d2w_dx_dt + d2w_dy_dt * mu))
    My = -D * (d2w_dx + mu * d2w_dx + tf * (d2w_dy_dt + d2w_dx_dt * mu))
    Mxy = D * (1 - mu) * (d2w_dx_dy + tf * d2w_dx_dy_dt)
    pprint("Mx")
    pprint(Mx)
    pprint("My")
    pprint(My)
    pprint("Mxy")
    pprint(Mxy)

tension_solve()
