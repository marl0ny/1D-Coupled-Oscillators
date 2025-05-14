"""Implementation of the recursive formula for obtaining the 
Hermite polynomials in an iterative manner.

This script will explicitly write out the first few of these polynomials
in C-style if or case statement code expressions.

Reference:

Shankar R., "The Harmonic Oscillator,"
in <i>Principles of Quantum Mechanics</i>, 2nd ed,
Springer, 1994, ch. 7., pg 185-221.

"""
from sympy import Symbol, simplify
import re


def factorial(n):
    res, prev = 1, 1
    for i in range(n+1):
        if i > 1:
            res = i*prev
            prev = res
    return res


def hermite(n, x):
    """Hermite polynomials.

    See Shankar, pg. 195, 7.3.21 to obtain the base cases,
    then 7.3.35 for the recursive relation itself.
    """
    prev1, prev2, res = 0, 0, 0
    for i in range(n+1):
        if i == 0:
            res = 1
            prev2 = 1
        elif i == 1:
            res = 2*x
            prev1 = 2*x
        elif i > 1:
            res = 2*(x*prev1 - (i-1)*prev2)
            prev2, prev1 = prev1, res
    return res


def get_number_and_end_position(expr, i):
    num = ''
    is_float = False
    while i < len(expr) and expr[i] in '0123456789eE.':
        if expr[i] == '.':
            is_float = True
            num += expr[i]
            i += 1
        if expr[i] in 'eE':
            is_float = True
            if (i + 1) < len(expr) and expr[i+1] in '+-':
                num += expr[i] + expr[i+1]
                i += 2
            else:
                num += expr[i]
                i += 1
        else:
            num += expr[i]
            i += 1
    if not is_float:
        num += '.0'
    return num, i


def replace_ints_to_floats(expr):
    new_expr = ''
    i = 0
    while (i < len(expr)):
        if expr[i] in '0123456789.' and (
            (i - 1 >= 0 and expr[i-1] in ' ()+-*/') 
            or i == 0):
            num, i = get_number_and_end_position(expr, i)
            new_expr += num
        else:
            new_expr += expr[i]
            i += 1
    return new_expr 


def write_case(i, x):
    expr = simplify(hermite(i, x))
    return f"case {str(i)}: " + "return "\
        + replace_ints_to_floats(
            re.sub("\*\*", '', str(expr)))\
        + ";"


def write_if_statement(i, x):
    expr = simplify(hermite(i, x))
    return f"{'else ' if i > 0 else ''}if (n == {str(i)}) " + "return "\
        + replace_ints_to_floats(
            re.sub("\*\*", '', str(expr)))\
        + ";"


if __name__ == '__main__':
    x_sym = Symbol('x')
    for i in range(20):
        print(write_if_statement(i, x_sym))

