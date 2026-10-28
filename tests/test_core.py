import unittest
from fractions import Fraction

from exact_polynomial_gcd_tools import Polynomial, polynomial_gcd


class TestPolynomial(unittest.TestCase):
    def test_zero_polynomial_degree(self):
        self.assertEqual(Polynomial([0]).degree(), -1)

    def test_zero_polynomial_is_zero(self):
        self.assertTrue(Polynomial([0]).is_zero())

    def test_trailing_zeros_stripped(self):
        self.assertEqual(Polynomial([1, 2, 0, 0]), Polynomial([1, 2]))

    def test_degree(self):
        self.assertEqual(Polynomial([1, 2, 3]).degree(), 2)

    def test_add(self):
        a = Polynomial([1, 2, 3])
        b = Polynomial([4, 5])
        self.assertEqual(a + b, Polynomial([5, 7, 3]))

    def test_sub(self):
        a = Polynomial([4, 5, 6])
        b = Polynomial([1, 2, 3])
        self.assertEqual(a - b, Polynomial([3, 3, 3]))

    def test_mul_by_scalar(self):
        a = Polynomial([1, 2, 3])
        self.assertEqual(a * 2, Polynomial([2, 4, 6]))
        self.assertEqual(2 * a, Polynomial([2, 4, 6]))

    def test_mul_by_fraction(self):
        a = Polynomial([1, 2, 3])
        self.assertEqual(a * Fraction(1, 2), Polynomial([Fraction(1, 2), 1, Fraction(3, 2)]))

    def test_mul_by_zero_scalar(self):
        a = Polynomial([1, 2, 3])
        self.assertEqual(a * 0, Polynomial([0]))

    def test_mul_by_polynomial(self):
        a = Polynomial([1, 1])
        b = Polynomial([1, 1])
        self.assertEqual(a * b, Polynomial([1, 2, 1]))

    def test_mul_by_zero_polynomial(self):
        a = Polynomial([1, 2, 3])
        self.assertEqual(a * Polynomial([0]), Polynomial([0]))

    def test_coeffs_returns_copy(self):
        a = Polynomial([1, 2, 3])
        c = a.coeffs
        c[0] = Fraction(99)
        self.assertEqual(a, Polynomial([1, 2, 3]))

    def test_equality(self):
        self.assertEqual(Polynomial([1, 2, 3]), Polynomial([1, 2, 3]))
        self.assertNotEqual(Polynomial([1, 2, 3]), Polynomial([1, 2]))

    def test_str_zero(self):
        self.assertEqual(str(Polynomial([0])), "0")

    def test_str_simple(self):
        self.assertEqual(str(Polynomial([2, 0, 3])), "3*x**2 + 2")

    def test_str_monic_linear(self):
        self.assertEqual(str(Polynomial([0, 1])), "x")


class TestGCD(unittest.TestCase):
    def test_coprime(self):
        a = Polynomial([1, 0, 1])  # x^2 + 1
        b = Polynomial([1, 1])      # x + 1
        self.assertEqual(polynomial_gcd(a, b), Polynomial([1]))

    def test_linear_factor(self):
        a = Polynomial([2, 1])                 # x + 2
        b = Polynomial([1, Fraction(1, 2)])    # (1/2)(x + 2)
        # Monic GCD is x + 2.
        self.assertEqual(polynomial_gcd(a, b), Polynomial([2, 1]))

    def test_quadratic_common_factor(self):
        # (x-1)(x-2) and (x-1)(x-3) share (x-1).
        a = Polynomial([2, -3, 1])  # (x-1)(x-2)
        b = Polynomial([3, -4, 1])  # (x-1)(x-3)
        self.assertEqual(polynomial_gcd(a, b), Polynomial([-1, 1]))

    def test_first_arg_zero(self):
        a = Polynomial([0])
        b = Polynomial([3, 2, 1])
        # GCD(0, b) = monic(b) = x^2 + 2x + 3.
        self.assertEqual(polynomial_gcd(a, b), Polynomial([3, 2, 1]))

    def test_second_arg_zero(self):
        a = Polynomial([3, 2, 1])
        b = Polynomial([0])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([3, 2, 1]))

    def test_both_zero(self):
        self.assertEqual(polynomial_gcd(Polynomial([0]), Polynomial([0])), Polynomial([0]))

    def test_scalar_gcd_is_one(self):
        a = Polynomial([5])
        b = Polynomial([7])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([1]))

    def test_non_integer_rationals(self):
        # a = (1/2)x + 1, b = (1/3)x + (2/3)
        # Both are scalar multiples of x + 2.
        a = Polynomial([1, Fraction(1, 2)])
        b = Polynomial([Fraction(2, 3), Fraction(1, 3)])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([2, 1]))

    def test_negative_leading_coefficient(self):
        # a = -(x^2 - 1) = -x^2 + 1, b = -(x + 1) = -x - 1
        # Monic GCD is x + 1.
        a = Polynomial([1, 0, -1])
        b = Polynomial([-1, -1])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([1, 1]))

    def test_gcd_is_monic(self):
        a = Polynomial([6, 3])  # 3(2 + x)
        b = Polynomial([4, 2])  # 2(2 + x)
        result = polynomial_gcd(a, b)
        self.assertEqual(result, Polynomial([2, 1]))
        self.assertEqual(result.coeffs[-1], Fraction(1))

    def test_type_error(self):
        with self.assertRaises(TypeError):
            polynomial_gcd("not a poly", Polynomial([1]))
        with self.assertRaises(TypeError):
            polynomial_gcd(Polynomial([1]), 42)

    def test_high_degree_coprime(self):
        # x^3 + 1 and x^3 - 1 are coprime over Q (their roots are disjoint).
        a = Polynomial([1, 0, 0, 1])
        b = Polynomial([-1, 0, 0, 1])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([1]))

    def test_repeated_factor(self):
        # a = (x-1)^2, b = (x-1)^3. GCD = (x-1)^2 = x^2 - 2x + 1.
        a = Polynomial([1, -2, 1])
        b = Polynomial([-1, 3, -3, 1])
        self.assertEqual(polynomial_gcd(a, b), Polynomial([1, -2, 1]))


if __name__ == "__main__":
    unittest.main()
