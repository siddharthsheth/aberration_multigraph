from math import comb, factorial
from numpy import prod

def count_for_c_l(n, c, l):
    # return prod([comb(n-1-sum(l[:i-1]), l[i]-1)*2**(l[i]-1)*factorial(l[i]-1) for i in range(c)])
    result = 1
    for i in range(c):
        result *= comb(n-1-sum(l[:i]), l[i]-1)*2**(l[i]-1)*factorial(l[i]-1)
    return result

def generate_l(n, c):
    m = n - 2*(c-1)
    l = [2]*(c-1).extend([m])
    yield l
    last = c-1
    while l[0] != m:
        while l[last] > 2:
            pass

# print(count_for_c_l(5, 2, (3,2)))