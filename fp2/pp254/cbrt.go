package pp254

import (
	"math/big"

	fp "github.com/yelhousni/fp2-cbrt/fields/pp254"
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
