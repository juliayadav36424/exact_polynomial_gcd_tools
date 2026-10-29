# exact_polynomial_gcd_tools

Computes the greatest common divisor of two polynomials with rational coefficients using the Euclidean algorithm over `fractions.Fraction`.

```python
from fractions import Fraction
from exact_polynomial_gcd_tools import Polynomial, polynomial_gcd

a = Polynomial([1, -3, 3, -1])      # (x-1)^3
b = Polynomial([1, -2, 1])          # (x-1)^2
print(polynomial_gcd(a, b))         # x**2 - 2*x + 1

# Rational coefficients are handled exactly:
c = Polynomial([Fraction(1, 2), 1])
d = Polynomial([Fraction(2, 3), Fraction(4, 3)])
print(polynomial_gcd(c, d))          # x + 2
```

## Why this exists

Symbolic algebra libraries exist and are excellent, but they bring a dependency graph the size of a small city. When the only thing you need is an exact GCD of two rational polynomials, `fractions.Fraction` plus the Euclidean algorithm is enough. The trade-off: no multivariate support, no modular arithmetic, no interpolation. In return you get a single file with zero third-party dependencies and a result you can trust bit-for-bit.

## The awkward edge

The GCD of two non-zero polynomials is returned in monic form (leading coefficient 1). The GCD of two zero polynomials is the zero polynomial itself, which has no leading coefficient to normalise; the function returns `Polynomial([0])` in that case. The degree of the zero polynomial is reported as -1, a convention chosen so the Euclidean loop terminates cleanly.
