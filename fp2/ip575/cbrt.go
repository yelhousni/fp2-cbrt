package ip575

import (
	"math/big"

	fp "github.com/yelhousni/fp2-cbrt/fields/ip575"
)

var lucasExponent = [9]uint64{
	12297829382473034411,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	12297829382473034410,
	1669334261878663850,
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

func (z *E2) cbrtVerifyAndAdjust(x *E2, y *E2) *E2 {
	var c E2
	c.Square(y).Mul(&c, y)
	if c.Equal(x) {
		return z.Set(y)
	}

	var omega = fp.Element{
		14017292452453342317,
		2409278055862628299,
		5947491282347698337,
		4737137772726116386,
		11574324834628588850,
		13362396080463313342,
		8790510260292158241,
		8998922594979982390,
		4866823027283296663,
	}
	var omega2 = fp.Element{
		4429451621256209294,
		16037466017846923316,
		12499252791361853278,
		13709606300983435229,
		6872419239080962765,
		5084347993246238273,
		9656233813417393374,
		9447821478729569225,
		1726446827187109480,
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
	for i := 571; i >= 1; i-- {
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
	var e0, e1 [9]uint64
	e0 = [9]uint64{16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 2225779015838218467}
	e1 = [9]uint64{4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 10248191152060862008, 4099276460824344803, 16397105843297379214, 556444753959554616}

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

	// 573 bits is odd: process top bit to align
	{
		b0 := (e0[8] >> 60) & 1
		b1 := (e1[8] >> 60) & 1
		idx := b0<<2 | b1
		if idx != 0 {
			z.Set(&table[idx-1])
		}
	}

	for i := 570; i >= 0; i -= 2 {
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
			if limb+1 < 9 {
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
