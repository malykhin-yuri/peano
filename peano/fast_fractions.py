from math import gcd
from functools import cached_property


class FastFraction:
    """Fast implementations of rational fractions, without implicit conversions."""

    # TODO: maybe cache fractions for -100<=n<=100, 1<=d<=100 ?

    def __init__(self, n, d):
        if d < 0:
            n = -n
            d = -d
        elif d == 0:
            raise ValueError("Division by zero!")
        g = gcd(n, d)
        self.n = n // g
        self.d = d // g

    @classmethod
    def convert(cls, x):
        if isinstance(x, cls):
            return x
        elif isinstance(x, int):
            return cls(x, 1)
        else:
            raise ValueError("Can't convert")

    @classmethod
    def parse(cls, fraction_str):
        if '/' in fraction_str:
            n, d = fraction_str.split('/')
        else:
            n, d = fraction_str, 1
        return cls(int(n), int(d))

    @cached_property
    def one_complement(self):
        return FastFraction(self.d - self.n, self.d)

    def __gt__(self, other):
        return self.n * other.d > other.n * self.d

    def __ge__(self, other):
        return self.n * other.d >= other.n * self.d

    def __lt__(self, other):
        return self.n * other.d < other.n * self.d

    def __le__(self, other):
        return self.n * other.d <= other.n * self.d

    def __eq__(self, other):
        return (self.n, self.d) == (other.n, other.d)

    def __neg__(self):
        return FastFraction(-self.n, self.d)

    def __mul__(self, other):
        return FastFraction(self.n * other.n, self.d * other.d)

    def __add__(self, other):
        return FastFraction(self.n * other.d + other.n * self.d, self.d * other.d)

    def __sub__(self, other):
        return FastFraction(self.n * other.d - other.n * self.d, self.d * other.d)

    def __truediv__(self, other):
        return FastFraction(self.n * other.d, self.d * other.n)

    def __pow__(self, power):
        return FastFraction(self.n**power, self.d**power)

    def __float__(self):
        return self.n / self.d

    # not compliant with float.int for negative values; gives lfloor
    def __int__(self):
        return self.n // self.d

    def __str__(self):
        if self.d == 1:
            return str(self.n)
        else:
            return '{}/{}'.format(self.n, self.d)

    def __hash__(self):
        return hash((self.n, self.d))

    def __repr__(self):
        return 'FastFraction({}, {})'.format(self.n, self.d)
