from __future__ import annotations

from fractions import Fraction
from typing import List


class Polynomial:
    """A polynomial in one variable with exact rational coefficients.

    Coefficients are stored in a list ``[c0, c1, c2, ...]`` where the term
    of degree ``i`` has coefficient ``c[i]``. The list is kept as short as
    possible: trailing zeros are stripped on construction so that two
    polynomials representing the same mathematical object always share the
    same internal representation. This makes equality and degree
    calculations trivial and keeps the Euclidean algorithm's termination
    condition easy to check.

    Using ``fractions.Fraction`` throughout means there is no rounding
    error at any step of the GCD computation. That is the whole point of
    this library: the result is exact, not an approximation.
    """

    __slots__ = ("_coeffs",)

    def __init__(self, coeffs):
        converted = [Fraction(c) for c in coeffs]
        # Strip trailing zeros so the representation is canonical.
        # We must keep a single zero when the polynomial is genuinely the
        # zero polynomial, otherwise degree() and is_zero() break.
        i = len(converted)
        while i > 1 and converted[i - 1] == 0:
            i -= 1
        self._coeffs = converted[:i]

    @property
    def coeffs(self) -> List[Fraction]:
        # Return a copy so callers cannot mutate our internal state.
        return list(self._coeffs)

    def degree(self) -> int:
        # Degree of the zero polynomial is defined as -1 here. This is a
        # convention; the Euclidean algorithm only needs to know that the
        # zero polynomial is "smaller" than every non-zero polynomial, and
        # -1 satisfies that.
        if len(self._coeffs) == 1 and self._coeffs[0] == 0:
            return -1
        return len(self._coeffs) - 1

    def is_zero(self) -> bool:
        return len(self._coeffs) == 1 and self._coeffs[0] == 0

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            return NotImplemented
        n = max(len(self._coeffs), len(other._coeffs))
        result = [Fraction(0)] * n
        for i, c in enumerate(self._coeffs):
            result[i] += c
        for i, c in enumerate(other._coeffs):
            result[i] += c
        return Polynomial(result)

    def __sub__(self, other):
        if not isinstance(other, Polynomial):
            return NotImplemented
        n = max(len(self._coeffs), len(other._coeffs))
        result = [Fraction(0)] * n
        for i, c in enumerate(self._coeffs):
            result[i] += c
        for i, c in enumerate(other._coeffs):
            result[i] -= c
        return Polynomial(result)

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            if self.is_zero() or other.is_zero():
                return Polynomial([0])
            result = [Fraction(0)] * (len(self._coeffs) + len(other._coeffs) - 1)
            for i, a in enumerate(self._coeffs):
                for j, b in enumerate(other._coeffs):
                    result[i + j] += a * b
            return Polynomial(result)
        if isinstance(other, (int, Fraction)):
            s = Fraction(other)
            if s == 0:
                return Polynomial([0])
            return Polynomial([c * s for c in self._coeffs])
        return NotImplemented

    __rmul__ = __mul__

    def __eq__(self, other) -> bool:
        if not isinstance(other, Polynomial):
            return NotImplemented
        return self._coeffs == other._coeffs

    def __hash__(self) -> int:
        return hash(tuple(self._coeffs))

    def __repr__(self) -> str:
        return f"Polynomial({self._coeffs!r})"

    def __str__(self) -> str:
        # Human-readable form. Not used by the algorithm, but helpful for
        # debugging and for error messages.
        if self.is_zero():
            return "0"
        parts = []
        for i in range(len(self._coeffs) - 1, -1, -1):
            c = self._coeffs[i]
            if c == 0:
                continue
            sign = "+" if c > 0 else "-"
            mag = abs(c)
            if i == 0:
                term = f"{mag}"
            elif i == 1:
                term = "x" if mag == 1 else f"{mag}*x"
            else:
                term = f"x**{i}" if mag == 1 else f"{mag}*x**{i}"
            parts.append((sign, term))
        sign, term = parts[0]
        out = term if sign == "+" else "-" + term
        for sign, term in parts[1:]:
            out += f" {sign} {term}"
        return out


def polynomial_gcd(a: Polynomial, b: Polynomial) -> Polynomial:
    """Return the greatest common divisor of two polynomials.

    The GCD is computed with the Euclidean algorithm over the rationals:
    repeatedly replace ``(a, b)`` with ``(b, a mod b)`` until ``b`` is the
    zero polynomial. The last non-zero remainder is the GCD.

    Because every coefficient is an exact ``Fraction``, no rounding ever
    occurs. The result is always returned with a leading coefficient of 1
    (monic normalisation); the only exception is the GCD of two zero
    polynomials, which is the zero polynomial itself and cannot be made
    monic.
    """
    if not isinstance(a, Polynomial):
        raise TypeError("a must be a Polynomial")
    if not isinstance(b, Polynomial):
        raise TypeError("b must be a Polynomial")

    x, y = a, b
    while not y.is_zero():
        x, y = y, _poly_mod(x, y)

    if x.is_zero():
        return Polynomial([0])

    # Normalise to monic form so the result is canonical: two callers
    # computing the GCD of the same pair always get the same answer.
    lead = x._coeffs[-1]
    if lead != 1:
        x = Polynomial([c / lead for c in x._coeffs])
    return x


def _poly_mod(a: Polynomial, b: Polynomial) -> Polynomial:
    """Return the remainder of ``a`` divided by ``b``.

    This is the workhorse of ``polynomial_gcd``. It performs polynomial
    long division. We do not need the quotient for the GCD algorithm, so
    this function discards it. ``b`` must be non-zero; the caller
    (``polynomial_gcd``) guarantees this via the ``while`` condition.
    """
    rem = a.coeffs
    bc = b._coeffs
    b_lead = bc[-1]
    b_deg = len(bc) - 1

    while True:
        if len(rem) < len(bc):
            break
        if len(rem) == 1 and rem[0] == 0:
            break
        # Cancel the leading term of the current remainder against the
        # leading term of b. This is one step of synthetic division.
        factor = rem[-1] / b_lead
        shift = len(rem) - len(bc)
        for i in range(len(bc)):
            rem[shift + i] -= factor * bc[i]
        # The leading coefficient is now zero by construction. Drop it so
        # the degree genuinely decreases each pass; otherwise the loop
        # would not terminate.
        rem.pop()

    if not rem:
        rem = [Fraction(0)]
    return Polynomial(rem)
