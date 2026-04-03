// Package fp2 contains per-prime implementations of Fp2 = Fp[u]/(u^2 - beta)
// arithmetic, including cube root computation via the algebraic torus T_2(Fp).
//
// # Torus cube root algorithm
//
// For x = x0 + x1*u in Fp2, the cube root is computed by reducing to base-field
// operations through the algebraic torus T_2(Fp):
//
//  1. Compute n = Norm(x) = x0^2 - beta*x1^2  and  m = cbrt(n) in Fp.
//  2. Map x to the torus via y = x^{q-1}, obtaining trace tau = Tr(y) = 2*(x0^2 - beta*x1^2)/n.
//  3. Cube-root on the torus using the Lucas V-sequence: sigma = V_e(tau) where e = 3^{-1} mod (q+1).
//  4. Recover the Fp2 result from sigma and m using Okeya-Sakurai-style formulae,
//     with Hamburg's trick to merge 1/U and 1/N into a single Fp exponentiation.
//
// # Edge cases: A1 = 0 and A0 = 0
//
// The general algorithm divides by x0*x1 and the norm, so it cannot handle the
// cases where one component is zero. These are treated separately. The question
// is whether falling back to a base-field cube root is complete, i.e., whether
// a non-cube in Fp could still have a cube root in Fp2.
//
// Expanding (a + b*u)^3 with u^2 = beta:
//
//	real part:  a^3 + 3*a*b^2*beta
//	imag part:  b*(3*a^2 + b^2*beta)
//
// For beta = -1 this simplifies to:
//
//	real:  a^3 - 3*a*b^2
//	imag:  b*(3*a^2 - b^2)
//
// A1 = 0 (x = x0, purely real):
//
// The imaginary part vanishes iff b = 0 or b^2 = 3*a^2.
//   - b = 0:        gives x0 = a^3, requiring x0 to be a cube in Fp.
//   - b^2 = 3*a^2:  gives x0 = a^3 - 9*a^3 = -8*a^3 = (-2a)^3,
//     which is always a perfect cube regardless of a.
//
// Therefore if x0 is not a cube in Fp, neither branch can produce it, and there
// is genuinely no cube root in Fp2. Returning nil is correct.
//
// A0 = 0 (x = x1*u, purely imaginary, beta = -1):
//
// The real part vanishes iff a = 0 or a^2 = 3*b^2.
//   - a = 0:        gives x1 = -b^3, requiring cbrt(-x1) in Fp.
//   - a^2 = 3*b^2:  gives x1 = 8*b^3, requiring cbrt(x1/8) in Fp.
//
// These two branches agree on cubicity because 8 = 2^3 is always a cube, so
// x1/8 is a cube iff x1 is, and -1 is always a cube mod any odd prime p with
// p = 1 mod 3 (since 6 | (p-1) implies (p-1)/3 is even, so (-1)^{(p-1)/3} = 1).
// When p = 2 mod 3 every element is a cube. Thus for all odd primes, -x1 is a
// cube iff x1 is a cube iff x1/8 is a cube, and both branches agree. The code
// only checks the a = 0 branch, which is sufficient.
//
// A0 = 0 with beta = -5 (e.g. BLS12-377 scalar field):
//
// The same reasoning applies: the a = 0 branch needs cbrt(-x1/5) and the
// a^2 = 15*b^2 branch needs cbrt(x1/40). Since -1 and 8 are always cubes,
// both reduce to whether x1/5 is a cube. The single base-field check is complete.
package fp2
