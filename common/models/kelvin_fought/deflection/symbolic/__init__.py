from sympy import *
from sympy import pi as PI

D = Symbol('D', real=True, nonzero=True)
pi = Symbol('rho_i', real=True, nonzero=True)
pw = Symbol('rho_w', real=True, nonzero=True)
h = Symbol('h', real=True, nonzero=True)
g = Symbol('g', real=True, nonzero=True)
t = Symbol('t', real=True, nonzero=True)
q = Symbol('q', real=True, nonzero=True)
r = Symbol('r', real=True, nonzero=True)
mu = Symbol('mu', real=True, nonzero=True)
H = Symbol('H', real=True, nonzero=True)
x = Symbol('x', real=True)
y = Symbol('y', real=True)
z = Symbol('z', real=True)
v = Symbol('v', real=True)
tf = Symbol('tau_phi', real=True, nonzero=True)
lm = Symbol('lambda', real=True, nonzero=True)
eta = Symbol('eta', real=True, nonzero=True)
delta = Symbol('delta', real=True)
PhiLE = Function('Phi')(lm, eta)
xi = Symbol('xi', real=True, nonzero=True)
chi = Symbol('chi', real=True, nonzero=True)
a = Symbol('a', real=True, nonzero=True)
b = Symbol('b', real=True, nonzero=True)
A = Symbol('A', real=True, nonzero=True)
B = Symbol('B', real=True, nonzero=True)
dPhidt = Symbol('dPhi_dt', real=True, nonzero=True)
P = Symbol('P', real=True, nonzero=True)
w1 = Function('w1')(x, y)
w2 = Function('w2')(t)
w = w1 + w2
Phi1 = Function('Phi1')(x, y)
Phi2 = Function('Phi2')(t)
Phi = Phi1 + Phi2


def nabla4(func):
    d4_dx4 = diff(func, x, 4)
    d4_dy4 = diff(func, y, 4)
    d4_dx2_dy2 = diff(diff(func, x, 2), y, 2)
    return d4_dx4 + 2 * d4_dx2_dy2 + d4_dy4


def diff_t(func):
    return diff(func, t) - v * diff(func, x)


def diff_t2(func):
    return diff(func, t, 2) - 2 * v * diff(diff(func, x), t) + v ** 2 * diff(func, x, 2)


delta = exp(-I * (lm * (x - v * t) + eta * y))


def fourier_integral(func):
    return 1 / (2 * PI) * Integral(func * delta, (lm, -oo, +oo), (eta, -oo, +oo))


model = D * (nabla4(w) + tf * diff_t(nabla4(w))) + pw * g * w + pi * h * diff_t2(w) + pw * diff_t(Phi) - P
model = simplify(model.subs(w2, 0).subs(Phi2, 0).doit())
pprint(model)


def deflection_solve(**specs):
    pass
    # if 'lm' in specs:
    #     lm = specs['lm']
    # else:
    #     lm = Symbol('lambda', real=True, nonzero=True)
    # if 'eta' in specs:
    #     eta = specs['eta']
    # else:
    #     eta = Symbol('eta', real=True, nonzero=True)
    # if lm == 0 and eta == 0:
    #     return 0, 0, 0
    # phi = lm * (v * t + xi - x) + eta * (chi - y)
    # delta = exp(I * phi)
    # Phi = Phile * cosh((H + z) * sqrt(lm ** 2 + eta ** 2)) * delta
    # dPhi_dz = diff(Phi, z)
    # wle = Symbol('omega_le')
    # w = wle * delta
    # border_eq = diff(w, t) - dPhi_dz.subs(z, 0)
    # Phi_le_solution = solve(border_eq, Phile)
    # try:
    #     Phi_le = Phi_le_solution[0]
    # except IndexError:
    #     return 0, 0, 0
    # print('hello')
    # Phi = simplify(Phi_le * cosh((H + z) * sqrt(lm ** 2 + eta ** 2)).subs(z, 0) * delta).doit()
    # dPhi_dt = diff(Phi, t)
    # model = D * (nabla4(w) + tf * diff(nabla4(w), t)) + pw * g * w + pi * h * diff(w, t,
    #                                                                                2) + pw * dPhidt + P * delta
    # wle = solve(model, wle)[0]
    #
    # K = tanh(H * sqrt(lm ** 2 + eta ** 2)) * sqrt(lm ** 2 + eta ** 2)
    # K.refine(Q.nonzero(K))
    # w = wle * delta
    # numer, denom = fraction(w)
    # numer = simplify(expand(P * numer, complex=True) / K)
    # denom = simplify(collect(simplify(denom), I).doit() / K)
    # A = collect(re(denom), D)
    # coefD = A.coeff(D)
    # A = (A - D * coefD) + D * factor(coefD)
    # B = factor(im(denom))
    #
    # dDelta_dxi = integrate(numer, xi, conds='none')
    # dDelta_dxi_dchi = integrate(dDelta_dxi, chi, conds='none')
    # numer = Subs(dDelta_dxi_dchi, xi, b).doit() - Subs(dDelta_dxi_dchi, xi, -b).doit()
    # numer = Subs(numer, chi, a).doit() - Subs(numer, chi, -a).doit()
    # Aa = Symbol('A', real=True, nonzero=True)
    # Bb = Symbol('B', real=True, nonzero=True)
    # conj = Aa - I * Bb
    # denom = simplify((Aa + I * Bb) * conj)
    # if lm == 0:
    #     numer = simplify(re(expand(numer * conj, complex=True)))
    #     return simplify(numer / denom), A, B
    # numer = - simplify(re(expand(numer * conj, complex=True)) * -1)
    # pprint(A)
    #
    # pprint(B)
    # return simplify(numer / denom), A, B
