package pp254

import (

	fp "github.com/yelhousni/fp2-cbrt/fp/pp254"
)

// lucasExponent is e = 3⁻¹ mod (p+1) as little-endian uint64 limbs.
var lucasExponent = [4]uint64{
	7593120314996402627,
	15936870763965662084,
	16724893366231265993,
	1162332755600990221,
}

// Cbrt sets z to the cube root of x and returns z.
// Returns nil if x is not a cubic residue.
func (z *E2) Cbrt(x *E2) *E2 {
	return z.cbrtTorus(x)
}

func (z *E2) cbrtTorus(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}

	if x.A0.IsZero() {
		var negA1 fp.Element
		negA1.Neg(&x.A1)
		var y E2
		if y.A1.Cbrt(&negA1) == nil {
			return nil
		}
		y.A0.SetZero()
		return z.cbrtVerifyAndAdjust(x, &y)
	}

	var x0sq, x1sq fp.Element
	x0sq.Square(&x.A0)
	x1sq.Square(&x.A1)
	// N = x₀² + x₁² (β = -1)
	var norm fp.Element
	norm.Add(&x0sq, &x1sq)

	m, normInv, ok := cbrtAndNormInverse(&norm)
	if !ok {
		return nil
	}

	// α_t = 2·(x₀² - β·x₁²)/N = 2·(x₀² + x₁²)/N ... wait, β=-1 so -β=1
	// Actually α_t = 2·(x₀² - β·x₁²)/N where β is the non-residue for u²=β
	// For BN254: u² = -1, so β = -1, and -β·x₁² = +x₁²
	// But that would give α_t = 2·(x₀² + x₁²)/N = 2·N/N = 2 always, which is wrong.
	//
	// The correct formula uses the trace of x^(p-1) on T₂:
	// If x = x₀ + x₁·u, then x^p = x₀ - x₁·u (Frobenius)
	// x^(p-1) = x^p / x = (x₀ - x₁u) / (x₀ + x₁u)
	// trace(x^(p-1)) = x^(p-1) + conj(x^(p-1)) = x^(p-1) + x^(1-p)
	//
	// For β = -1: N(x) = x₀² + x₁², trace = 2(x₀² - x₁²)/N
	// (note: -β = -(-1) = +1 → x₀² + β·x₁² → no, let's re-derive)
	//
	// trace(x^(p-1)) = 2·Re(x^p · x̄) / |x|² = 2·(x₀² - β·x₁²) / N
	// With β = -1: 2·(x₀² - (-1)·x₁²)/N = 2·(x₀² + x₁²)/N = 2
	// This means α_t = 2 always when β = -1 and x₀,x₁ both nonzero?
	// No! The Re(x^p · x̄) = x₀² + x₁² only when Frobenius acts as conjugation.
	// Let me re-check:
	// x^p = x₀ + x₁·u^p. Since u² = -1, u^p = u^(p mod 2).
	// p is odd, so u^p = u·(u²)^((p-1)/2) = u·(-1)^((p-1)/2)
	// For BN254, (p-1)/2 is odd (p ≡ 3 mod 4), so u^p = -u
	// Thus x^p = x₀ - x₁·u = conjugate(x)
	// trace(x^(p-1)) = x^(p-1) + x^(-(p-1)) = x^(p-1) + (x^(p-1))^(-1) · since on T₂
	// Wait, x^(p-1) has norm 1, so its inverse IS its conjugate.
	// Actually on T₂: if r is on the torus, r·r̄ = 1.
	// r = x^(p-1) = x^p/x = x̄/x
	// trace(r) = r + r̄ = x̄/x + x/x̄
	// = (x̄² + x²)/(x·x̄) = (x̄² + x²)/N(x)
	// x̄ = x₀ - x₁u, x = x₀ + x₁u
	// x̄² = x₀² + x₁²β - 2x₀x₁u = (x₀² - x₁²) - 2x₀x₁u  [β=-1]
	// x² = (x₀² - x₁²) + 2x₀x₁u
	// x̄² + x² = 2(x₀² - x₁²)
	// So trace = 2(x₀² - x₁²)/N, where N = x₀² + x₁²
	// This is NOT always 2! My mistake above.

	var alphaT fp.Element
	alphaT.Sub(&x0sq, &x1sq) // x₀² - x₁²
	alphaT.Double(&alphaT)
	alphaT.Mul(&alphaT, &normInv)

	sp := lucasV(&alphaT)

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

// cbrtVerifyAndAdjust verifies y³ = x, adjusting by 9th roots of unity if needed.
// For p mod 9 = 1: there are primitive 9th roots of unity in Fp.
func (z *E2) cbrtVerifyAndAdjust(x *E2, y *E2) *E2 {
	var c E2
	c.Square(y).Mul(&c, y)
	if c.Equal(x) {
		return z.Set(y)
	}

	var omega = fp.Element{
		8183898218631979349,
		12014359695528440611,
		12263358156045030468,
		3187210487005268291,
	}
	var omega2 = fp.Element{
		3697675806616062876,
		9065277094688085689,
		6918009208039626314,
		2775033306905974752,
	}
	var zeta = fp.Element{
		9092840637269024442,
		11284133545212953584,
		7919372827184455520,
		1596114425137527684,
	}
	var zeta2 = fp.Element{
		1735008219140503419,
		10465829585049341007,
		6017168831245289042,
		1570250484855163800,
	}

	var cw2 E2
	cw2.A0.Mul(&c.A0, &omega2)
	cw2.A1.Mul(&c.A1, &omega2)
	if cw2.Equal(x) {
		z.A0.Mul(&y.A0, &zeta)
		z.A1.Mul(&y.A1, &zeta)
		return z
	}

	var cw E2
	cw.A0.Mul(&c.A0, &omega)
	cw.A1.Mul(&c.A1, &omega)
	if cw.Equal(x) {
		z.A0.Mul(&y.A0, &zeta2)
		z.A1.Mul(&y.A1, &zeta2)
		return z
	}

	return nil
}

// cbrtAndNormInverse for BN254: p mod 9 = 1
// t = norm^((q-19)/27), m = norm · t, normInv = m^17 · t^10
func cbrtAndNormInverse(norm *fp.Element) (m, normInv fp.Element, ok bool) {
	var t fp.Element
	t.ExpByCbrtHelperQMinus19Div27(*norm)
	m.Mul(norm, &t)

	// normInv = m^17 · t^10
	var m2 fp.Element
	m2.Square(&m)
	var t2, t4, t8 fp.Element
	t2.Square(&t)
	t4.Square(&t2)
	t8.Square(&t4)
	normInv.Mul(&t8, &t2) // t^10

	var m4, m8, m16 fp.Element
	m4.Square(&m2)
	m8.Square(&m4)
	m16.Square(&m8)
	t2.Mul(&m16, &m) // m^17 (reuse t2)
	normInv.Mul(&normInv, &t2)

	// Verify m³ = norm, adjust by ζ if needed
	var c fp.Element
	c.Mul(&m2, &m)
	if !c.Equal(norm) {
		var zeta = fp.Element{
			9092840637269024442,
			11284133545212953584,
			7919372827184455520,
			1596114425137527684,
		}
		var zeta2 = fp.Element{
			1735008219140503419,
			10465829585049341007,
			6017168831245289042,
			1570250484855163800,
		}
		var omega = fp.Element{
			8183898218631979349,
			12014359695528440611,
			12263358156045030468,
			3187210487005268291,
		}
		var omega2 = fp.Element{
			3697675806616062876,
			9065277094688085689,
			6918009208039626314,
			2775033306905974752,
		}

		var cw2 fp.Element
		cw2.Mul(&c, &omega2)
		if cw2.Equal(norm) {
			m.Mul(&m, &zeta)
		} else {
			var cw fp.Element
			cw.Mul(&c, &omega)
			if cw.Equal(norm) {
				m.Mul(&m, &zeta2)
			} else {
				return m, normInv, false
			}
		}
	}

	return m, normInv, true
}

func lucasV(alpha *fp.Element) fp.Element {
	var v0, v1, two fp.Element
	two.SetUint64(2)
	v0.Set(alpha)
	v1.Square(alpha).Sub(&v1, &two)

	var prod fp.Element
	for i := 252 - 1; i >= 1; i-- {
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
	v0.Mul(&v0, &v1).Sub(&v0, alpha)
	return v0
}

// CbrtDirect computes the cube root via direct Fp2 exponentiation.
func (z *E2) cbrtDirect(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}
	var y E2
	y.expByE2Cbrt(*x)
	return z.cbrtVerifyAndAdjust(x, &y)
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
	var y E2
	y.expByE2CbrtFrobenius(*x)
	return z.cbrtVerifyAndAdjust(x, &y)
}

func (z *E2) expByE2CbrtFrobenius(x E2) *E2 {
	var e0, e1 [4]uint64
	e0 = [4]uint64{11330118615407775612, 13329482909165958675, 6094068723619673846, 1420628923512321382}
	e1 = [4]uint64{3736998300411372985, 15839356218909848207, 7815919431097959468, 258296167911331160}

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

	// 253 bits is odd: process top bit to align
	{
		b0 := (e0[3] >> 60) & 1
		b1 := (e1[3] >> 60) & 1
		idx := b0<<2 | b1
		if idx != 0 {
			z.Set(&table[idx-1])
		}
	}

	for i := 250; i >= 0; i -= 2 {
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
			if limb+1 < 4 {
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
		var negA1 fp.Element
		negA1.Neg(&x.A1)
		var y E2
		if y.A1.Cbrt(&negA1) == nil {
			return nil
		}
		y.A0.SetZero()
		return z.cbrtVerifyAndAdjust(x, &y)
	}

	var x0sq, x1sq fp.Element
	x0sq.Square(&x.A0)
	x1sq.Square(&x.A1)
	var norm fp.Element
	norm.Add(&x0sq, &x1sq)

	m, normInv, ok := cbrtAndNormInverse(&norm)
	if !ok {
		return nil
	}

	var alphaT fp.Element
	alphaT.Sub(&x0sq, &x1sq)
	alphaT.Double(&alphaT)
	alphaT.Mul(&alphaT, &normInv)

	sp := lucasVPrac(&alphaT)

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

func (z *E2) cbrtFrobenius1bit(x *E2) *E2 {
	if x.A1.IsZero() {
		if z.A0.Cbrt(&x.A0) == nil {
			return nil
		}
		z.A1.SetZero()
		return z
	}
	var y E2
	y.expByE2CbrtFrobenius1bit(*x)
	return z.cbrtVerifyAndAdjust(x, &y)
}

func (z *E2) expByE2CbrtFrobenius1bit(x E2) *E2 {
	var e0, e1 [4]uint64
	e0 = [4]uint64{11330118615407775612, 13329482909165958675, 6094068723619673846, 1420628923512321382}
	e1 = [4]uint64{3736998300411372985, 15839356218909848207, 7815919431097959468, 258296167911331160}

	var xFrob E2
	xFrob.Conjugate(&x)

	var table [3]E2
	table[0].Set(&xFrob)      // (0,1)
	table[1].Set(&x)          // (1,0)
	table[2].Mul(&x, &xFrob)  // (1,1)

	z.SetOne()

	for i := 253 - 1; i >= 0; i-- {
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

// expByE2Cbrt is equivalent to z.Exp(x, 0xad76de41a5af60ea315d97688a2286dda209e0149bcd7bb708258216df4faba1830243f61450cc77484dacc5bf5165ad7b68d3edc5898e503f230287a81acb).
// It uses an addition chain generated by github.com/mmcloughlin/addchain.
//
// Operations: 497 squares 99 multiplies
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
		t18 E2
		t19 E2
		t20 E2
		t21 E2
		t22 E2
		t23 E2
		t24 E2
		t25 E2
		t26 E2
		t27 E2
		t28 E2
		t29 E2
		t30 E2
		t31 E2
		t32 E2
		t33 E2
		t34 E2
		t35 E2
		t36 E2
		t37 E2
		t38 E2
		t39 E2
		t40 E2
		t41 E2
	)

	t0.Square(&x)
	t1.Square(&t0)
	t2.Mul(&x, &t1)
	t3.Mul(&x, &t2)
	t4.Mul(&t1, &t3)
	t5.Mul(&x, &t4)
	t6.Mul(&t3, &t5)
	t7.Mul(&x, &t6)
	t8.Mul(&t1, &t6)
	t9.Mul(&t5, &t7)
	t10.Mul(&t3, &t9)
	t11.Mul(&t3, &t10)
	t12.Mul(&t1, &t11)
	t13.Mul(&t0, &t12)
	t14.Mul(&t7, &t11)
	t15.Mul(&t1, &t14)
	t16.Mul(&t1, &t15)
	t17.Mul(&t1, &t16)
	t18.Mul(&t1, &t17)
	t19.Mul(&t1, &t18)
	t20.Mul(&t0, &t19)
	t21.Mul(&t4, &t19)
	t22.Mul(&t7, &t21)
	t23.Mul(&t1, &t22)
	t24.Mul(&t3, &t23)
	t25.Mul(&t7, &t23)
	t26.Mul(&t0, &t25)
	t27.Mul(&t0, &t26)
	t28.Mul(&t0, &t27)
	t29.Mul(&t0, &t28)
	t30.Mul(&t7, &t29)
	t31.Mul(&t3, &t30)
	t3.Mul(&t7, &t31)
	t32.Mul(&t1, &t3)
	t33.Mul(&t1, &t32)
	t34.Mul(&t0, &t33)
	t35.Mul(&t1, &t34)
	t36.Mul(&t1, &t35)
	t37.Mul(&t0, &t36)
	t38.Mul(&t4, &t37)
	t39.Mul(&t0, &t38)
	t0.Mul(&t1, &t39)
	t40.Mul(&t4, &t0)
	t4.Mul(&t7, &t0)
	t7.Mul(&t1, &t4)
	t1.Mul(&t26, &t0)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Mul(&t7, &t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Square(&t41)
	t41.Square(&t1)
	t1.Mul(&t13, &t41)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t26, &t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Mul(&t18, &t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Mul(&t34, &t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t26, &t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t8, &t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t36, &t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t14, &t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Mul(&t33, &t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Square(&t1)
	t1.Square(&t13)
	t13.Mul(&t6, &t1)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t20, &t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t31, &t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t32, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t38, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t19, &t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t11, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t23, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t0, &t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t33, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Mul(&t40, &t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t18, &t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Square(&t13)
	t13.Square(&t6)
	t6.Mul(&t27, &t13)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Mul(&t32, &t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Mul(&t39, &t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Mul(&t24, &t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Mul(&t9, &t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Mul(&t35, &t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Square(&t6)
	t6.Square(&t27)
	t27.Mul(&t25, &t6)
	t25.Square(&t27)
	t27.Square(&t25)
	t25.Square(&t27)
	t27.Square(&t25)
	t25.Square(&t27)
	t27.Square(&t25)
	t25.Square(&t27)
	t27.Square(&t25)
	t25.Square(&t27)
	t27.Square(&t25)
	t25.Mul(&t28, &t27)
	t28.Square(&t25)
	t25.Square(&t28)
	t28.Square(&t25)
	t25.Square(&t28)
	t28.Square(&t25)
	t25.Square(&t28)
	t28.Mul(&t14, &t25)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Mul(&t20, &t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Mul(&t16, &t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Square(&t14)
	t14.Square(&t28)
	t28.Mul(&t37, &t14)
	t37.Square(&t28)
	t28.Square(&t37)
	t37.Square(&t28)
	t28.Square(&t37)
	t37.Square(&t28)
	t28.Square(&t37)
	t37.Square(&t28)
	t28.Square(&t37)
	t37.Square(&t28)
	t28.Mul(&t4, &t37)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Mul(&t30, &t28)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Mul(&t3, &t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Mul(&t32, &t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Mul(&t24, &t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Square(&t4)
	t4.Square(&t30)
	t30.Mul(&t21, &t4)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Mul(&t22, &t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Mul(&t34, &t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Square(&t21)
	t21.Square(&t30)
	t30.Mul(&t12, &t21)
	t12.Square(&t30)
	t30.Square(&t12)
	t12.Square(&t30)
	t30.Square(&t12)
	t12.Square(&t30)
	t30.Square(&t12)
	t12.Square(&t30)
	t30.Square(&t12)
	t12.Square(&t30)
	t30.Square(&t12)
	t12.Square(&t30)
	t30.Mul(&t39, &t12)
	t39.Square(&t30)
	t30.Square(&t39)
	t39.Square(&t30)
	t30.Square(&t39)
	t39.Square(&t30)
	t30.Square(&t39)
	t39.Square(&t30)
	t30.Square(&t39)
	t39.Mul(&t7, &t30)
	t7.Square(&t39)
	t39.Square(&t7)
	t7.Square(&t39)
	t39.Square(&t7)
	t7.Square(&t39)
	t39.Square(&t7)
	t7.Square(&t39)
	t39.Square(&t7)
	t7.Mul(&t36, &t39)
	t36.Square(&t7)
	t7.Square(&t36)
	t36.Square(&t7)
	t7.Square(&t36)
	t36.Square(&t7)
	t7.Square(&t36)
	t36.Square(&t7)
	t7.Square(&t36)
	t36.Mul(&t29, &t7)
	t29.Square(&t36)
	t36.Square(&t29)
	t29.Square(&t36)
	t36.Square(&t29)
	t29.Square(&t36)
	t36.Square(&t29)
	t29.Square(&t36)
	t36.Mul(&t17, &t29)
	t29.Square(&t36)
	t36.Square(&t29)
	t29.Square(&t36)
	t36.Square(&t29)
	t29.Square(&t36)
	t36.Mul(&t2, &t29)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Mul(&t15, &t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Mul(&t10, &t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Mul(&t31, &t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Mul(&t24, &t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Square(&t36)
	t36.Square(&t2)
	t2.Mul(&t22, &t36)
	t22.Square(&t2)
	t2.Square(&t22)
	t22.Square(&t2)
	t2.Square(&t22)
	t22.Square(&t2)
	t2.Square(&t22)
	z.Mul(&t5, &t2)

	return z
}
