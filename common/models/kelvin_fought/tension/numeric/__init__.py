import multiprocessing as mp

import numpy as np

from common.models.kelvin_fought.deflection.numeric import deflection_func as w, integrator_adapter


# def f1(a):
#     return - a['l'] ** 2 * w(a)
#
#
# def f2(a):
#     return - a['e'] ** 2 * w(a)


def f3(a):
    return a['l'] * a['e'] * w(a)


# def f4(a):
#     return - a['l'] ** 3 * a['v'] * w(a)


# def f5(a):
#     return


def f6(a):
    return - a['l'] ** 2 * a['v'] * w(a)


def dmx(a):
    w_val = w(a)
    f1 = - a['l'] ** 2 * w_val
    f2 = - a['e'] ** 2 * w_val
    f4 = - a['l'] ** 3 * a['v'] * w_val
    f5 = - - a['e'] ** 2 * a['v'] * w_val
    return - a['D'] * (f1 + a['mu'] * f2 + a['t_f'] * (f4(a) + a['mu'] * f5(a)))


def dmy(a):
    return - a['D'] * (f2(a) + a['mu'] * f1(a) + a['t_f'] * (f5(a) + a['mu'] * f4(a)))


def dmxy(a):
    return - a['D'] * (1 - a['mu']) * (f3(a) + a['t_f'] * f6(a))


def tension_func_sx(a):
    val = 6 * dmx(a) / a['h'] ** 2
    return val


def tension_func_sy(a):
    return 6 * dmy(a) / a['h'] ** 2


def tension_func_txy(a):
    return 6 * dmxy(a) / a['h'] ** 2


def tension_state(a, xrange, yrange, specs):
    pool = mp.Pool(processes=mp.cpu_count())
    data_sx = None
    data_sy = None
    data_txy = None
    if 'sx' in specs:
        results = dict()
        data_sx = dict()
        for y in yrange:
            results[y] = [pool.apply_async(func=integrate_for, args=(x, y, tension_func_sx, a.copy())) for x in xrange]
        for y in results:
            xlist = results[y]
            for x in xlist:
                x_p, y_p, v = x.get()
                data_sx[x_p, y_p] = v
    if 'sy' in specs:
        results = dict()
        data_sy = dict()
        for y in yrange:
            results[y] = [pool.apply_async(func=integrate_for, args=(x, y, tension_func_sy, a.copy())) for x in xrange]
        for y in results:
            xlist = results[y]
            for x in xlist:
                x_p, y_p, v = x.get()
                data_sy[x_p, y_p] = v
    if 'txy' in specs:
        results = dict()
        data_txy = dict()
        for y in yrange:
            results[y] = [pool.apply_async(func=integrate_for, args=(x, y, tension_func_txy, a.copy())) for x in xrange]
        for y in results:
            xlist = results[y]
            for x in xlist:
                x_p, y_p, v = x.get()
                data_txy[x_p, y_p] = v
    return data_sx, data_sy, data_txy


def integrate_for(x, y, func, a):
    a['x'] = x
    a['y'] = y
    a['v'] *= 0.99
    a['t'] = 0
    print(str(x) + ';' + str(y))
    return x, y, -4 * a['P'] * integrator_adapter(func, a, 0, np.inf, 0, np.inf) / 1000 / (np.pi ** 2)
