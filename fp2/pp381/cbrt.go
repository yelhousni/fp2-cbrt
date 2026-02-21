package pp381

import (
	"math/big"

	fp "github.com/yelhousni/fp2-cbrt/fields/pp381"
)

var lucasExponent = [6]uint64{
	10616391696595805071,
	736713837172402858,
	14776387573661061644,
	14710942035944605247,
	1804034592823567431,
	624599539215846622,
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
		14772873186050699377,
		6749526151121446354,
		6372666795664677781,
		10283423008382700446,
		286397964926079186,
		1796971870900422465,
	}
	var omega2 = fp.Element{
		3526659474838938856,
		17562030475567847978,
		1632777218702014455,
		14009062335050482331,
		3906511377122991214,
		368068849512964448,
	}
	var zeta = fp.Element{
		13616190144799058984,
		9227582506135211912,
		4426607408274926740,
		7455198167498346307,
		10794825842164118204,
		335101026345095675,
	}
	var zeta2 = fp.Element{
		3828863564860874189,
		5918733612565202776,
		16843310164143221096,
		16127847466718491017,
		17435063908385505950,
		407112797415018074,
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

// cbrtAndNormInverse for BLS12-381: p mod 9 = 1
// t = norm^((q-10)/27), m = norm · t², normInv = m^8 · t^11
func cbrtAndNormInverse(norm *fp.Element) (m, normInv fp.Element, ok bool) {
	var t fp.Element
	t.ExpByCbrtHelperQMinus10Div27(*norm)
	var t2 fp.Element
	t2.Square(&t)
	m.Mul(norm, &t2)

	var m2 fp.Element
	m2.Square(&m)
	var m4, m8 fp.Element
	m4.Square(&m2)
	m8.Square(&m4)

	var t4, t8 fp.Element
	t4.Square(&t2)
	t8.Square(&t4)
	normInv.Mul(&t8, &t2)
	normInv.Mul(&normInv, &t) // t^11
	normInv.Mul(&normInv, &m8)

	var c fp.Element
	c.Mul(&m2, &m)
	if !c.Equal(norm) {
		var zeta = fp.Element{
			13616190144799058984,
			9227582506135211912,
			4426607408274926740,
			7455198167498346307,
			10794825842164118204,
			335101026345095675,
		}
		var zeta2 = fp.Element{
			3828863564860874189,
			5918733612565202776,
			16843310164143221096,
			16127847466718491017,
			17435063908385505950,
			407112797415018074,
		}
		var omega = fp.Element{
			14772873186050699377,
			6749526151121446354,
			6372666795664677781,
			10283423008382700446,
			286397964926079186,
			1796971870900422465,
		}
		var omega2 = fp.Element{
			3526659474838938856,
			17562030475567847978,
			1632777218702014455,
			14009062335050482331,
			3906511377122991214,
			368068849512964448,
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
	for i := 379 - 1; i >= 1; i-- {
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
	var e0, e1 [6]uint64
	e0 = [6]uint64{5647076082758821762, 4917847391015903535, 16418208415178957382, 12246214690225216582, 16351950493800281736, 693999488017607357}
	e1 = [6]uint64{13477428459872568307, 4181133553843500676, 1641820841517895738, 15982016727990162951, 14547915900976714304, 69399948801760735}

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

	for i := 378; i >= 0; i -= 2 {
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
