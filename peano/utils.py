"""
Ratio and utility functions.
"""

from collections import Counter
from functools import lru_cache
import itertools
import math

from quicktions import Fraction


def ratio_linf(d, dv, dt):
    return Fraction(max(abs(x) for x in dv) ** d, dt)


def ratio_l1(d, dv, dt):
    return Fraction(sum(abs(x) for x in dv) ** d, dt)


def ratio_l2_squared(d, dv, dt):
    return Fraction(sum(x ** 2 for x in dv) ** d, dt ** 2)


def ratio_l2(d, dv, dt):
    assert d % 2 == 0
    d2 = d // 2
    return Fraction(sum(x ** 2 for x in dv) ** d2, dt)


@lru_cache(maxsize=2**20)
def get_int_cube_with_cache(dim, div, cubes):
    """Integer coordinates for sequence of embedded cubes."""
    x = [0] * dim
    div_power = 1
    for cube in reversed(cubes):
        for j in range(dim):
            x[j] += cube[j] * div_power
        div_power *= div
    return x


@lru_cache(maxsize=2**20)
def get_int_time_with_cache(dim, div, cnums):
    """Integer time for sequence of cnums of embedded cubes."""
    G = div**dim
    # multiply by G**l, l = len(cnums), i.e. depth
    # t = c0/G + c1/G**2 + ... = (c_{l-1} + c_{l-2}*G + ..) / G^l
    t = 0
    Gpower = 1
    for cnum in reversed(cnums):
        t += cnum * Gpower
        Gpower *= G
    return t


def get_periodic_sum(start, period, d):
    """
    Sum the non-periodic and periodic parts.

    The sum is:
    s_0/d + s_1/d^2 + ... + s_{k-1}/d^k (non-periodic part = start) +
       + s_k/d^{k+1} + ... + s_{k+m-1}/d^{k+m} + ... (periodic part = period)
    """
    n0 = 0
    d_power = 1
    for x in reversed(start):
        n0 += x * d_power
        d_power *= d
    t0 = Fraction(n0, d_power)

    np = 0
    dp_power = 1
    for x in reversed(period):
        np += x * dp_power
        dp_power *= d
    tp = Fraction(np, d_power * (dp_power - 1))

    return t0 + tp


def get_lcm(iterable):
    """Least common multiple of integer sequence."""
    lcm = 1
    for x in iterable:
        lcm = (lcm * x) // math.gcd(lcm, x)
    return lcm


def combinations_product(iter_ids, iter_dict):
    """
    Product of combinations of iterables.
    Args:
        iter_ids  --  list of ids
        iter_dict --  dict {iter_id: iterable}

    Generate tuples (elem_1,...,elem_n) where elem_j comes from iter_ids[j] iterable.
    All positions with given iter_id are equivalent, so
    we take combinations of iterables within that group, not whole product.
    """
    counts = Counter(iter_ids)
    combs = [itertools.combinations_with_replacement(iter_dict[iter_id], r=cnt) for iter_id, cnt in counts.items()]
    id2idx = {iter_id: idx for idx, iter_id in enumerate(counts)}
    for groups_items in itertools.product(*combs):
        groups_items = [list(elems) for elems in groups_items]
        result = []
        for iter_id in iter_ids:
            result.append(groups_items[id2idx[iter_id]].pop(0))
        yield tuple(result)
