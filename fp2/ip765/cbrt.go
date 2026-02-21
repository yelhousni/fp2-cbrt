package ip765

import (
	"math/big"

	fp "github.com/yelhousni/fp2-cbrt/fields/ip765"
)

// lucasExponent is e = 3⁻¹ mod (p+1) as little-endian uint64 limbs.
var lucasExponent = [12]uint64{
	12297829382473034411,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	385808368078072490,
}

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

// cbrtVerifyAndAdjust for p mod 9 = 4:
// 3 | (p-1) so cube roots of unity exist, but 9 ∤ (p-1) so no 9th roots.
func (z *E2) cbrtVerifyAndAdjust(x *E2, y *E2) *E2 {
	var c E2
	c.Square(y).Mul(&c, y)
	if c.Equal(x) {
		return z.Set(y)
	}

	// Primitive cube roots of unity ω, ω² in Fp (Montgomery form)
	var omega = fp.Element{
		5038930390546469911,
		17787860233557347803,
		9849330006372172977,
		14474249416463688568,
		15913942290039186735,
		189484831851572713,
		8634920383567828977,
		16880187712022066618,
		11181262368701633890,
		8765771123670726457,
		3891436893190388257,
		490302666823792294,
	}
	var omega2 = fp.Element{
		13407813683163081688,
		658883840152203812,
		8597414067337378638,
		3972494657245863047,
		2532801783670364880,
		18257259241857978902,
		9811823690141722638,
		1566556361687484997,
		7265481705007917725,
		9680972950038825158,
		14555307180519163358,
		739180031448353113,
	}

	var cw2 E2
	cw2.A0.Mul(&c.A0, &omega2)
	cw2.A1.Mul(&c.A1, &omega2)
	if cw2.Equal(x) {
		z.A0.Mul(&y.A0, &omega)
		z.A1.Mul(&y.A1, &omega)
		var check E2
		check.Square(z).Mul(&check, z)
		if check.Equal(x) {
			return z
		}
	}

	var cw E2
	cw.A0.Mul(&c.A0, &omega)
	cw.A1.Mul(&c.A1, &omega)
	if cw.Equal(x) {
		z.A0.Mul(&y.A0, &omega2)
		z.A1.Mul(&y.A1, &omega2)
		var check E2
		check.Square(z).Mul(&check, z)
		if check.Equal(x) {
			return z
		}
	}

	return nil
}

// cbrtAndNormInverse for p mod 9 = 4:
// t = norm^((q-4)/9), m = norm · t², normInv = m² · t⁵
func cbrtAndNormInverse(norm *fp.Element) (m, normInv fp.Element, ok bool) {
	var t fp.Element
	t.ExpByCbrtHelperQMinus4Div9(*norm)
	var t2 fp.Element
	t2.Square(&t)
	m.Mul(norm, &t2) // m = norm · t²

	// normInv = m² · t⁵
	var m2 fp.Element
	m2.Square(&m)

	var t4 fp.Element
	t4.Square(&t2)
	normInv.Mul(&t4, &t) // t⁵
	normInv.Mul(&normInv, &m2)

	// Verify m³ = norm
	var c fp.Element
	c.Mul(&m2, &m)
	if !c.Equal(norm) {
		return m, normInv, false
	}

	return m, normInv, true
}

func lucasV(alpha *fp.Element) fp.Element {
	var v0, v1, two fp.Element
	two.SetUint64(2)
	v0.Set(alpha)
	v1.Square(alpha).Sub(&v1, &two)

	var prod fp.Element
	for i := 761; i >= 1; i-- {
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
	var e0, e1 [12]uint64
	e0 = [12]uint64{16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 514411157437429987}
	e1 = [12]uint64{4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 128602789359357496}

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

	// 763 bits is odd: process top bit to align
	{
		b0 := (e0[11] >> 58) & 1
		b1 := (e1[11] >> 58) & 1
		idx := b0<<2 | b1
		if idx != 0 {
			z.Set(&table[idx-1])
		}
	}

	for i := 760; i >= 0; i -= 2 {
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
			if limb+1 < 12 {
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
