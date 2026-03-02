package pp377

import (

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
		if z.A1.Cbrt(&negA1OverBeta) == nil {
			return nil
		}
		z.A0.SetZero()
		return z.cbrtVerify(x)
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

	// τ = 2·(x₀² - |β|·x₁²)/N = trace of x^{p-1} on T₂
	var tau fp.Element
	tau.Sub(&x0sq, &betaX1sq)
	tau.Double(&tau)
	tau.Mul(&tau, &normInv)

	// σ = V_e(τ) where e = 3⁻¹ mod (p+1)
	sigma := lucasV(&tau)

	// Recovery: z₀ = x₀/(m·(σ-1)), z₁ = x₁/(m·(σ+1))
	var one, d0, d1, d0d1, d0d1Inv fp.Element
	one.SetOne()
	d0.Sub(&sigma, &one)
	d0.Mul(&m, &d0)
	d1.Add(&sigma, &one)
	d1.Mul(&m, &d1)

	d0d1.Mul(&d0, &d1)
	d0d1Inv.Inverse(&d0d1)

	z.A0.Mul(&d1, &d0d1Inv).Mul(&z.A0, &x.A0)
	z.A1.Mul(&d0, &d0d1Inv).Mul(&z.A1, &x.A1)

	return z.cbrtVerify(x)
}

// cbrtVerify checks z³ = x. Returns nil if not.
// For p mod 9 = 7, cube root is unique: no adjustment needed.
func (z *E2) cbrtVerify(x *E2) *E2 {
	var c E2
	c.Square(z).Mul(&c, z)
	if !c.Equal(x) {
		return nil
	}
	return z
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
func (z *E2) cbrtDirect(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}

	z.expByE2Cbrt(*x)
	return z.cbrtVerify(x)
}

// CbrtFrobenius computes the E2 cube root using Frobenius decomposition
// and 2-bit windowed Strauss-Shamir multi-exponentiation.
func (z *E2) cbrtFrobenius(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}
	z.expByE2CbrtFrobenius(*x)
	return z.cbrtVerify(x)
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

func (z *E2) cbrtTorusPrac(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}

	if x.A0.IsZero() {
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
		if z.A1.Cbrt(&negA1OverBeta) == nil {
			return nil
		}
		z.A0.SetZero()
		return z.cbrtVerify(x)
	}

	var x0sq, x1sq fp.Element
	x0sq.Square(&x.A0)
	x1sq.Square(&x.A1)
	var norm, betaX1sq fp.Element
	betaX1sq.Set(&x1sq)
	fp.MulBy5(&betaX1sq)
	norm.Add(&x0sq, &betaX1sq)

	m, normInv := cbrtAndNormInverse(&norm)

	var tau fp.Element
	tau.Sub(&x0sq, &betaX1sq)
	tau.Double(&tau)
	tau.Mul(&tau, &normInv)

	sigma := lucasVPrac(&tau)

	var one, d0, d1, d0d1, d0d1Inv fp.Element
	one.SetOne()
	d0.Sub(&sigma, &one)
	d0.Mul(&m, &d0)
	d1.Add(&sigma, &one)
	d1.Mul(&m, &d1)

	d0d1.Mul(&d0, &d1)
	d0d1Inv.Inverse(&d0d1)

	z.A0.Mul(&d1, &d0d1Inv).Mul(&z.A0, &x.A0)
	z.A1.Mul(&d0, &d0d1Inv).Mul(&z.A1, &x.A1)

	return z.cbrtVerify(x)
}

func (z *E2) cbrtFrobenius1bit(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}
	z.expByE2CbrtFrobenius1bit(*x)
	return z.cbrtVerify(x)
}

func (z *E2) expByE2CbrtFrobenius1bit(x E2) *E2 {
	var e0, e1 [6]uint64
	e0 = [6]uint64{17623453223078942038, 922513019478125226, 15586486611553363968, 3095920025918078116, 16134117903398705782, 67276840392497054}
	e1 = [6]uint64{14428078918715397461, 369005207791250090, 2545245829879435264, 8617065639851051893, 17521693605585213282, 26910736156998821}

	var xFrob E2
	xFrob.Conjugate(&x)

	var table [3]E2
	table[0].Set(&xFrob)      // (0,1)
	table[1].Set(&x)          // (1,0)
	table[2].Mul(&x, &xFrob)  // (1,1)

	z.SetOne()

	for i := 376 - 1; i >= 0; i-- {
		z.Square(z)

		limb := i / 64
		bit := uint(i % 64)
		b0 := (e0[limb] >> bit) & 1
		b1 := (e1[limb] >> bit) & 1

		idx := b0<<1 | b1
		if idx != 0 {
			z.Mul(z, &table[idx-1])
		}
	}

	return z
}

// expByE2Cbrt is equivalent to z.Exp(x, 0xa0ac6746271b5cf6caf0d2875dd4829ad3aa2498a59c3fdc192af14cf4bf1ea0057c62016ad86392de0411d082fb9dc97b4257efc987277279f04c87e65bdf10e2a07dc3f1e965f796ca7c063000199ad968355555559075aaaaaaaaaaab).
// It uses an addition chain generated by github.com/mmcloughlin/addchain.
//
// Operations: 747 squares 138 multiplies
func (z *E2) expByE2Cbrt(x E2) *E2 {
	var (
		t0 E2
		t1 E2
		t2 E2
		t3 E2
		t4 E2
		t5 E2
		t6 E2
		t7 E2
		t8 E2
		t9 E2
		t10 E2
		t11 E2
		t12 E2
		t13 E2
		t14 E2
		t15 E2
		t16 E2
		t17 E2
	)

	t0.Square(&x)
	t1.Mul(&x, &t0)
	t2.Mul(&t0, &t1)
	t3.Mul(&t0, &t2)
	t4.Mul(&t0, &t3)
	t5.Mul(&t0, &t4)
	t6.Mul(&t0, &t5)
	t7.Mul(&t0, &t6)
	t8.Mul(&t0, &t7)
	t9.Mul(&t0, &t8)
	t10.Mul(&t0, &t9)
	t11.Mul(&t0, &t10)
	t12.Mul(&t0, &t11)
	t13.Mul(&t0, &t12)
	t14.Mul(&t0, &t13)
	t15.Mul(&t0, &t14)
	t0.Square(&t15)
	t16.Mul(&x, &t0)
	t0.Mul(&t8, &t16)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t10, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t8, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t9, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t2, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t1, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t9, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t8, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t5, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t11, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t7, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t13, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t10, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t3, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t6, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t2, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t14, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t14, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t6, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t4, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t2, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t6, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t6, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t14, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t2, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t4, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t9, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t2, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t5, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t3, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t16, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t1, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t3, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t12, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t10, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t7, &t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t2, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t12, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t14, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&x, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t16, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Mul(&t7, &t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Square(&t17)
	t17.Square(&t0)
	t0.Mul(&t2, &t17)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t10, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t7, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t1, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&x, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t5, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t10, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t5, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t1, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t3, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t4, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t13, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t1, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&x, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t8, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t6, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&x, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t11, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t13, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t9, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t11, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t4, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t7, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t6, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Mul(&t4, &t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&x, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t16, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t16, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t9, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t3, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t9, &t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Square(&t2)
	t2.Square(&t0)
	t0.Mul(&t11, &t2)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&t9, &t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&t12, &t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Mul(&t7, &t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&t9, &t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&x, &t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Mul(&t16, &t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Mul(&t12, &t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&t13, &t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Square(&t11)
	t11.Square(&t0)
	t0.Mul(&t13, &t11)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Mul(&t3, &t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Mul(&x, &t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Mul(&t3, &t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Mul(&t10, &t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Mul(&t15, &t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Mul(&t3, &t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Square(&t13)
	t13.Square(&t0)
	t0.Mul(&t16, &t13)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Mul(&t7, &t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Mul(&t4, &t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Mul(&t12, &t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Mul(&t15, &t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Mul(&t7, &t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Mul(&t5, &t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Mul(&t12, &t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Mul(&t9, &t0)
	t0.Square(&t16)
	t16.Square(&t0)
	t0.Square(&t16)
	t16.Mul(&t3, &t0)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Mul(&t1, &t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Mul(&t1, &t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Mul(&t12, &t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Mul(&t9, &t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Square(&t16)
	t16.Square(&t3)
	t3.Mul(&t5, &t16)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t12, &t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t6, &t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t6, &t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Mul(&t10, &t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t10, &t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Mul(&t10, &t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t10, &t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Mul(&t12, &t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Square(&t3)
	t3.Square(&t5)
	t5.Mul(&t14, &t3)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Mul(&t6, &t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Mul(&t10, &t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Mul(&t10, &t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Mul(&t10, &t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Mul(&t10, &t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Mul(&t10, &t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Mul(&t10, &t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Square(&t5)
	t5.Square(&t14)
	t14.Mul(&t10, &t5)
	t10.Square(&t14)
	t14.Square(&t10)
	t10.Square(&t14)
	z.Mul(&t1, &t10)

	return z
}
