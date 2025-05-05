from sympy import Symbol, simplify

def factorial(n):
    res, prev = 1, 1
    for i in range(n+1):
        if i > 1:
            res = i*prev
            prev = res
    return res


def hermite(n, x):
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

x = Symbol('x')
for i in range(10):
    print(str(simplify(hermite(i, x))))

# for i in range(6):
#     print(factorial(i))

