package pp381

import (

	fp "github.com/yelhousni/fp2-cbrt/fp/pp381"
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

	var tau fp.Element
	tau.Sub(&x0sq, &x1sq)
	tau.Double(&tau)
	tau.Mul(&tau, &normInv)

	sigma := lucasV(&tau)

	var one, sigmaM1, sigmaP1, d0, d1, d0d1, d0d1Inv fp.Element
	one.SetOne()
	sigmaM1.Sub(&sigma, &one)
	sigmaP1.Add(&sigma, &one)
	d0.Mul(&m, &sigmaM1)
	d1.Mul(&m, &sigmaP1)

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

	var tau fp.Element
	tau.Sub(&x0sq, &x1sq)
	tau.Double(&tau)
	tau.Mul(&tau, &normInv)

	sigma := lucasVPrac(&tau)

	var one, sigmaM1, sigmaP1, d0, d1, d0d1, d0d1Inv fp.Element
	one.SetOne()
	sigmaM1.Sub(&sigma, &one)
	sigmaP1.Add(&sigma, &one)
	d0.Mul(&m, &sigmaM1)
	d1.Mul(&m, &sigmaP1)

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
	var e0, e1 [6]uint64
	e0 = [6]uint64{5647076082758821762, 4917847391015903535, 16418208415178957382, 12246214690225216582, 16351950493800281736, 693999488017607357}
	e1 = [6]uint64{13477428459872568307, 4181133553843500676, 1641820841517895738, 15982016727990162951, 14547915900976714304, 69399948801760735}

	var xFrob E2
	xFrob.Conjugate(&x)

	var table [3]E2
	table[0].Set(&xFrob)      // (0,1)
	table[1].Set(&x)          // (1,0)
	table[2].Mul(&x, &xFrob)  // (1,1)

	z.SetOne()

	for i := 380 - 1; i >= 0; i-- {
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

// expByE2Cbrt is equivalent to z.Exp(x, 0x190b8ad76f8849c0701770fc867ca9d8feb0087bcb44fd3337e96b01f2e8bbdd0fa2d9f75d8c3cff998773ab047aa139fa626e17edf07656dbcc0fb8513ed34fa847c66a9bea57d169eef1e7300bbd895e206963317cfcdb818e38e49be8d3).
// It uses an addition chain generated by github.com/mmcloughlin/addchain.
//
// Operations: 751 squares 141 multiplies
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
	)

	t0.Square(&x)
	t1.Mul(&x, &t0)
	t2.Mul(&x, &t1)
	t3.Mul(&x, &t2)
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
	t16.Mul(&t0, &t15)
	t17.Mul(&t0, &t16)
	t18.Mul(&t0, &t17)
	t19.Mul(&t2, &t18)
	t2.Mul(&t0, &t19)
	t20.Mul(&t0, &t2)
	t21.Mul(&t0, &t20)
	t22.Mul(&t0, &t21)
	t23.Mul(&t0, &t22)
	t24.Mul(&t0, &t23)
	t25.Mul(&t0, &t24)
	t26.Mul(&t0, &t25)
	t27.Mul(&t0, &t26)
	t28.Mul(&t0, &t27)
	t29.Mul(&t0, &t28)
	t30.Mul(&t0, &t29)
	t0.Mul(&t19, &t29)
	t31.Mul(&t14, &t0)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Mul(&t12, &t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Mul(&t20, &t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Mul(&t12, &t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Mul(&t26, &t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Mul(&t23, &t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Square(&t0)
	t0.Square(&t32)
	t32.Mul(&t5, &t0)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t4, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t4, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t12, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t4, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t30, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t17, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t19, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t13, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t2, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t14, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t31, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t6, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t17, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t28, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t13, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t7, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t19, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t15, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t24, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t26, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t15, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t21, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t1, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t16, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t12, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t9, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t28, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t26, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t17, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t29, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t21, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t19, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t26, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t12, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t23, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t17, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t4, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t31, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t1, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t24, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t28, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t19, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t20, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t18, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t25, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t17, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t27, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t16, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t10, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t10, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t4, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t22, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Mul(&t14, &t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t16, &t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Square(&t32)
	t32.Square(&t5)
	t5.Mul(&t28, &t32)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t20, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t14, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t8, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t1, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t16, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t4, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t3, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t19, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t14, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t10, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t19, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Mul(&t25, &t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t18, &t5)
	t5.Square(&t28)
	t28.Square(&t5)
	t5.Square(&t28)
	t28.Mul(&t4, &t5)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t24, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Mul(&t11, &t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t26, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Mul(&t25, &t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t20, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t15, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t21, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Mul(&t29, &t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t26, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Mul(&t18, &t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t24, &t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Square(&t4)
	t4.Square(&t28)
	t28.Mul(&t10, &t4)
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
	t4.Square(&t28)
	t28.Mul(&t12, &t4)
	t12.Square(&t28)
	t28.Square(&t12)
	t12.Square(&t28)
	t28.Square(&t12)
	t12.Square(&t28)
	t28.Square(&t12)
	t12.Square(&t28)
	t28.Mul(&t29, &t12)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Mul(&t9, &t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Mul(&t20, &t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Mul(&t23, &t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Mul(&t7, &t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Square(&t29)
	t29.Square(&t28)
	t28.Mul(&t6, &t29)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Mul(&t24, &t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Mul(&t22, &t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Square(&t28)
	t28.Square(&t6)
	t6.Mul(&t19, &t28)
	t19.Square(&t6)
	t6.Square(&t19)
	t19.Square(&t6)
	t6.Square(&t19)
	t19.Square(&t6)
	t6.Square(&t19)
	t19.Mul(&t27, &t6)
	t27.Square(&t19)
	t19.Square(&t27)
	t27.Square(&t19)
	t19.Square(&t27)
	t27.Square(&t19)
	t19.Square(&t27)
	t27.Mul(&t21, &t19)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Mul(&t1, &t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Mul(&t23, &t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Mul(&t23, &t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Mul(&t23, &t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Square(&t21)
	t21.Square(&t27)
	t27.Mul(&t13, &t21)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Mul(&t10, &t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Mul(&t16, &t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Mul(&t18, &t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	t27.Square(&t13)
	t13.Square(&t27)
	z.Mul(&t10, &t13)

	return z
}
