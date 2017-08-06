import multiprocessing as mp
from builtins import print
from sys import float_info

import numpy as np
from scipy.integrate import nquad

eps = float_info[8]


def deflection_func(a):
    k = np.sqrt(a['l'] ** 2 + a['e'] ** 2)
    A = a['D'] * k ** 4 + a['g'] * a['p_w'] - 3 * a['h'] * a['l'] ** 2 * a['p_i'] * a['v'] ** 2
    A -= 4 * (a['l'] ** 2 * a['v'] ** 2 * a['p_w']) / (k * np.tanh(a['H'] * k))
    B = a['D'] * a['l'] * a['t_f'] * a['v'] * (k ** 4 + a['e'] ** 4)
    phi = a['e'] * a['y'] + a['l'] * (a['x'] - a['v'] * a['t'])
    C = np.sin(a['e'] * a['a']) * np.sin(a['l'] * a['b'])
    return (A * np.cos(phi) - B * np.sin(phi)) * C / ((A ** 2 + B ** 2) * a['l'] * a['e'])
    # # a['v'] = 0
    # if a['l'] == 0 and a['e'] == 0:
    #     return 0
    # l_sq = a['l'] ** 2
    # e_sq = a['e'] ** 2
    # le_sq = l_sq + e_sq
    # le_sqrt = sqrt(le_sq)
    # if a['l'] == 0:
    #     numer = a['b'] * sin(a['a'] * a['e']) * cos(a['y'] * a['e'])
    #     denom = a['e'] * (a['D'] * e_sq ** 2 + a['g'] * a['p_w'])
    #     return numer / denom
    # if a['e'] == 0:
    #     A = a['D'] * l_sq ** 2 + a['g'] * a['p_w'] - a['h'] * l_sq * a['p_i'] * a['v'] ** 2 - (
    #         (l_sq * a['p_w'] * a['v'] ** 2) / tanh(a['H'] * fabs(a['l'])) * fabs(a['l']))
    #     B = a['D'] * a['l'] ** 5 * a['t_f'] * a['v']
    #     phi = a['l'] * (a['x'] - a['v'] * a['t'])
    #     numer = a['a'] * (A * cos(phi) - B * sin(phi)) * sin(a['b'] * a['l'])
    #     denom = a['l'] * (A ** 2 + B ** 2)
    #     return numer / denom
    # A = a['D'] * le_sq ** 2 + a['p_w'] * a['g'] - l_sq * a['v'] ** 2 * (a['p_i'] * a['h'] + a['p_w'] / (
    #     tanh(a['H'] * le_sqrt) * le_sqrt))
    # B = a['D'] * a['t_f'] * a['l'] * a['v'] * le_sq ** 2
    # phi = a['e'] * a['y'] + a['l'] * (a['x'] - a['v'] * a['t'])
    # numer = (A * cos(phi) - B * sin(phi)) * sin(a['e'] * a['a']) * sin(a['l'] * a['b'])
    # denom = a['l'] * a['e'] * (A ** 2 + B ** 2)
    # return numer / denom


def deflection_state(func, a, xrange, yrange):
    pool = mp.Pool(processes=mp.cpu_count())
    results = dict()
    data = dict()
    for y in yrange:
        results[y] = [pool.apply_async(func=integrate_for, args=(x, y, func, a.copy())) for x in xrange]
    for y in results:
        xlist = results[y]
        for x in xlist:
            x_p, y_p, v = x.get()
            data[x_p, y_p] = v
    return data


def integrator_adapter(func, args, inner_l, inner_u, outer_l, outer_u):
    def R2Func(x, y):
        args['l'] = x
        args['e'] = y
        return func(args)

    return nquad(lambda x, y: R2Func(x, y), [[inner_l, inner_u], [outer_l, outer_u]], full_output=True)[0]


def integrate_for(x, y, func, a):
    a['x'] = x
    a['y'] = y
    a['v'] -= 0.1
    print(str(x) + ';' + str(y))
    # print(a)
    return x, y, -16 * a['P'] * integrator_adapter(func, a, 0, np.inf, 0, np.inf) / (np.pi * 2)
