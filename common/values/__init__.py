from numpy import sqrt


class Values:
    __a = {
        't_f': 10,
        'p_i': 900.0,
        'g': 9.81,
        'h': 0.5,
        'E': 5e9,
        'mu': 0.3,
        'H': 5.0,
        'p_w': 1000.0,
        'a': 15,
        'b': 6.5,
        'P': 400000,
        'x': 0,
        'y': 0,
        'l': 0,
        'e': 0,
        't': 0,
        'v': None,
        'D': None
    }

    __forced = False

    def __init__(self, forced=False):
        self.__forced = forced
        a = self.__a
        a['v'] = sqrt(a['g'] * a['H'])
        a['D'] = a['E'] * a['h'] ** 3 / (12 * (1 - a['mu'] ** 2))

    def __getitem__(self, item):
        return self.__a[item]

    def __setitem__(self, key, value):
        forced = self.__forced
        a = self.__a
        if key == 'D':
            if forced:
                a[key] = value
            else:
                pass
        elif key == 'E' or key == 'h' or key == 'mu':
            a[key] = value
            if not forced:
                a['D'] = a['E'] * a['h'] ** 3 / (12 * (1 - a['mu'] ** 2))
        elif key == 'g' or key == 'H':
            a[key] = value
            if not forced:
                a['v'] = sqrt(a['g'] * a['H'])
        else:
            a[key] = value

    def copy(self):
        c = Values()
        c.__a = self.__a.copy()
        c.__forced = self.__forced
        return c

    def __str__(self):
        return self.__a.__str__() + ' (f = ' + str(self.__forced) + ')'

    def set_forced(self, forced):
        self.__forced = forced

    def get_forced(self):
        return self.__forced
