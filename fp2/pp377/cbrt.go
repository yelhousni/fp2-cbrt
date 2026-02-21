package pp377

import (
	"math/big"

	fp "github.com/yelhousni/fp2-cbrt/fp/pp377"
)

// lucasExponent is e = 3⁻¹ mod (p+1) as little-endian uint64 limbs.
var lucasExponent = [6]uint64{
	3195374304363544577,
	553507811686875136,
	13041240781673928704,
	12925598459776577839,
	17059168371523044115,
	40366104235498232,
}

// Cbrt sets z to the cube root of x and returns z.
// Returns nil if x is not a cubic residue.
// For BLS12-377, p mod 9 = 7, so the cube root is unique.
func (z *E2) Cbrt(x *E2) *E2 {
	return z.cbrtTorus(x)
}

// cbrtTorus computes the cube root of x in E2 using the algebraic torus T₂(Fp).
func (z *E2) cbrtTorus(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}

	if x.A0.IsZero() {
		// Fp2 = Fp[u]/(u² - (-5)), so u³ = -5·u
		// x = x₁·u → y = a·u → y³ = a³·(-5)·u → a³ = x₁/(-5)
		var negA1OverBeta fp.Element
		betaInvNeg := fp.Element{
			330620507644336508,
			9878087358076053079,
			11461392860540703536,
			6973035786057818995,
			8846909097162646007,
			104838758629667239,
		}
		negA1OverBeta.Neg(&x.A1)
		negA1OverBeta.Mul(&negA1OverBeta, &betaInvNeg)
		var y E2
		if y.A1.Cbrt(&negA1OverBeta) == nil {
			return nil
		}
		y.A0.SetZero()
		return z.cbrtVerifyAndAdjust(x, &y)
	}

	// x₀², x₁²
	var x0sq, x1sq fp.Element
	x0sq.Square(&x.A0)
	x1sq.Square(&x.A1)
	// N = x₀² + 5·x₁² (norm of x, since beta = -5)
	var norm, betaX1sq fp.Element
	betaX1sq.Set(&x1sq)
	fp.MulBy5(&betaX1sq)
	norm.Add(&x0sq, &betaX1sq)

	// m = cbrt(N) and normInv = 1/N
	m, normInv := cbrtAndNormInverse(&norm)

	// α_t = 2·(x₀² - β·x₁²)/N = trace of x^{p-1} on T₂
	var alphaT fp.Element
	alphaT.Sub(&x0sq, &betaX1sq)
	alphaT.Double(&alphaT)
	alphaT.Mul(&alphaT, &normInv)

	// s₁ = V_e(α_t, 1) where e = 3⁻¹ mod (p+1)
	sp := lucasV(&alphaT)

	// Recovery: z₀ = x₀/(m·(s₁-1)), z₁ = x₁/(m·(s₁+1))
	var one, s1m1, s1p1, d0, d1, d0d1, d0d1Inv fp.Element
	one.SetOne()
	s1m1.Sub(&sp, &one)
	s1p1.Add(&sp, &one)
	d0.Mul(&m, &s1m1)
	d1.Mul(&m, &s1p1)

	d0d1.Mul(&d0, &d1)
	d0d1Inv.Inverse(&d0d1)

	var y E2
	y.A0.Mul(&d1, &d0d1Inv).Mul(&y.A0, &x.A0)
	y.A1.Mul(&d0, &d0d1Inv).Mul(&y.A1, &x.A1)

	return z.cbrtVerifyAndAdjust(x, &y)
}

// cbrtVerifyAndAdjust verifies y³ = x. Returns nil if not.
// For p mod 9 = 7, cube root is unique: no adjustment needed.
func (z *E2) cbrtVerifyAndAdjust(x *E2, y *E2) *E2 {
	var c E2
	c.Square(y).Mul(&c, y)
	if !c.Equal(x) {
		return nil
	}
	return z.Set(y)
}

// cbrtAndNormInverse computes m = cbrt(norm) and normInv = 1/norm.
// For p mod 9 = 7: t = norm^((q-7)/9), m = norm · t, normInv = m⁵ · t⁴
func cbrtAndNormInverse(norm *fp.Element) (m, normInv fp.Element) {
	var t fp.Element
	t.ExpByCbrtHelperQMinus7Div9(*norm)

	// m = norm · t
	m.Mul(norm, &t)

	// normInv = m^5 · t^4
	var m2, m4 fp.Element
	m2.Square(&m)
	m4.Square(&m2)
	var mPow fp.Element
	mPow.Mul(&m4, &m) // m⁵

	var t2, t4 fp.Element
	t2.Square(&t)
	t4.Square(&t2)

	normInv.Mul(&mPow, &t4)

	return m, normInv
}

// lucasV computes V_e(alpha, 1) where e = 3⁻¹ mod (p+1).
func lucasV(alpha *fp.Element) fp.Element {
	var v0, v1, two fp.Element
	two.SetUint64(2)
	v0.Set(alpha)
	v1.Square(alpha).Sub(&v1, &two)

	var prod fp.Element

	for i := 375 - 1; i >= 1; i-- {
		bit := (lucasExponent[i/64] >> uint(i%64)) & 1
		prod.Mul(&v0, &v1).Sub(&prod, alpha)
		if bit == 0 {
			v1.Set(&prod)
			v0.Square(&v0).Sub(&v0, &two)
		} else {
			v0.Set(&prod)
			v1.Square(&v1).Sub(&v1, &two)
		}
	}

	// Last bit (bit 0) is 1
	v0.Mul(&v0, &v1).Sub(&v0, alpha)
	return v0
}

// CbrtDirect computes the cube root using direct E2 exponentiation for benchmarking comparison.
// Uses e = 3⁻¹ mod ((p²-1)/3) as the exponent.
func (z *E2) CbrtDirect(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}

	p := fp.Modulus()
	n := new(big.Int).Mul(p, p)
	n.Sub(n, big.NewInt(1)) // p²-1
	// Remove all factors of 3 from n
	three := big.NewInt(3)
	for new(big.Int).Mod(n, three).Sign() == 0 {
		n.Div(n, three)
	}
	e := new(big.Int).ModInverse(three, n)

	var y E2
	y.Exp(*x, e)
	return z.cbrtVerifyAndAdjust(x, &y)
}

// CbrtFrobenius computes the E2 cube root using Frobenius decomposition
// and 2-bit windowed Strauss-Shamir multi-exponentiation.
func (z *E2) CbrtFrobenius(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}
	var y E2
	y.expByE2CbrtFrobenius(*x)
	return z.cbrtVerifyAndAdjust(x, &y)
}

func (z *E2) expByE2CbrtFrobenius(x E2) *E2 {
	var e0, e1 [6]uint64
	e0 = [6]uint64{17623453223078942038, 922513019478125226, 15586486611553363968, 3095920025918078116, 16134117903398705782, 67276840392497054}
	e1 = [6]uint64{14428078918715397461, 369005207791250090, 2545245829879435264, 8617065639851051893, 17521693605585213282, 26910736156998821}

	var xFrob E2
	xFrob.Conjugate(&x)

	var table [15]E2
	var xFrob2, xFrob3, x2, x3 E2
	x2.Square(&x)
	x3.Mul(&x2, &x)
	xFrob2.Square(&xFrob)
	xFrob3.Mul(&xFrob2, &xFrob)

	table[0].Set(&xFrob)
	table[1].Set(&xFrob2)
	table[2].Set(&xFrob3)
	table[3].Set(&x)
	table[4].Mul(&x, &xFrob)
	table[5].Mul(&x, &xFrob2)
	table[6].Mul(&x, &xFrob3)
	table[7].Set(&x2)
	table[8].Mul(&x2, &xFrob)
	table[9].Mul(&x2, &xFrob2)
	table[10].Mul(&x2, &xFrob3)
	table[11].Set(&x3)
	table[12].Mul(&x3, &xFrob)
	table[13].Mul(&x3, &xFrob2)
	table[14].Mul(&x3, &xFrob3)

	z.SetOne()

	for i := 374; i >= 0; i -= 2 {
		z.Square(z)
		z.Square(z)

		limb := i / 64
		bit := uint(i % 64)

		var w0, w1 uint64
		if bit < 63 {
			w0 = (e0[limb] >> bit) & 3
			w1 = (e1[limb] >> bit) & 3
		} else {
			w0 = (e0[limb] >> 63) & 1
			w1 = (e1[limb] >> 63) & 1
			if limb+1 < 6 {
				w0 |= (e0[limb+1] & 1) << 1
				w1 |= (e1[limb+1] & 1) << 1
			}
		}

		idx := w0<<2 | w1
		if idx != 0 {
			z.Mul(z, &table[idx-1])
		}
	}

	return z
}
