// Fp2 = Fp[u]/(u² + 5) arithmetic for BLS12-377
package pp377

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/fp2-cbrt/fp/pp377"
)

// E2 is a degree-two extension of fp.Element
type E2 struct {
	A0, A1 fp.Element
}

func (z *E2) SetZero() *E2   { z.A0.SetZero(); z.A1.SetZero(); return z }
func (z *E2) SetOne() *E2    { z.A0.SetOne(); z.A1.SetZero(); return z }
func (z *E2) Set(x *E2) *E2  { z.A0.Set(&x.A0); z.A1.Set(&x.A1); return z }
func (z *E2) Equal(x *E2) bool { return z.A0.Equal(&x.A0) && z.A1.Equal(&x.A1) }
func (z *E2) IsZero() bool   { return z.A0.IsZero() && z.A1.IsZero() }

func (z *E2) SetRandom() (*E2, error) {
	if _, err := z.A0.SetRandom(); err != nil {
		return nil, err
	}
	if _, err := z.A1.SetRandom(); err != nil {
		return nil, err
	}
	return z, nil
}

func (z *E2) Add(x, y *E2) *E2 {
	z.A0.Add(&x.A0, &y.A0)
	z.A1.Add(&x.A1, &y.A1)
	return z
}

func (z *E2) Sub(x, y *E2) *E2 {
	z.A0.Sub(&x.A0, &y.A0)
	z.A1.Sub(&x.A1, &y.A1)
	return z
}

func (z *E2) Double(x *E2) *E2 {
	z.A0.Double(&x.A0)
	z.A1.Double(&x.A1)
	return z
}

func (z *E2) Neg(x *E2) *E2 {
	z.A0.Neg(&x.A0)
	z.A1.Neg(&x.A1)
	return z
}

func (z *E2) Conjugate(x *E2) *E2 {
	z.A0.Set(&x.A0)
	z.A1.Neg(&x.A1)
	return z
}

// Mul sets z = x * y in Fp2 = Fp[u]/(u² + 5)
// Uses Karatsuba: (a+bu)(c+du) = (ac - 5bd) + (ad+bc)u
func (z *E2) Mul(x, y *E2) *E2 {
	var a, b, c fp.Element
	a.Add(&x.A0, &x.A1)
	b.Add(&y.A0, &y.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &y.A0)
	c.Mul(&x.A1, &y.A1)
	z.A1.Sub(&a, &b).Sub(&z.A1, &c)
	fp.MulBy5(&c)
	z.A0.Sub(&b, &c)
	return z
}

// Square sets z = x² in Fp2 = Fp[u]/(u² + 5)
func (z *E2) Square(x *E2) *E2 {
	var c0, c2 fp.Element
	c0.Add(&x.A0, &x.A1)
	c2.Neg(&x.A1)
	fp.MulBy5(&c2)
	c2.Add(&c2, &x.A0)

	c0.Mul(&c0, &c2)
	c2.Mul(&x.A0, &x.A1).Double(&c2)
	z.A1 = c2
	c2.Double(&c2)
	z.A0.Add(&c0, &c2)

	return z
}

// Inverse sets z = 1/x in Fp2
func (z *E2) Inverse(x *E2) *E2 {
	var t0, t1, tmp fp.Element
	a := &x.A0
	b := &x.A1
	t0.Square(a)
	t1.Square(b)
	tmp.Set(&t1)
	fp.MulBy5(&tmp)
	t0.Add(&t0, &tmp)
	t1.Inverse(&t0)
	z.A0.Mul(a, &t1)
	z.A1.Mul(b, &t1).Neg(&z.A1)
	return z
}

// norm sets x to the norm of z: N(z) = z₀² + 5·z₁²
func (z *E2) norm(x *fp.Element) {
	var tmp fp.Element
	x.Square(&z.A1)
	tmp.Set(x)
	fp.MulBy5(&tmp)
	x.Square(&z.A0).Add(x, &tmp)
}

// MulByElement multiplies z by a scalar in Fp
func (z *E2) MulByElement(x *E2, y *fp.Element) *E2 {
	var yCopy fp.Element
	yCopy.Set(y)
	z.A0.Mul(&x.A0, &yCopy)
	z.A1.Mul(&x.A1, &yCopy)
	return z
}

// Exp sets z = x^k (mod q²)
func (z *E2) Exp(x E2, k *big.Int) *E2 {
	if k.IsUint64() && k.Uint64() == 0 {
		return z.SetOne()
	}
	e := k
	if k.Sign() == -1 {
		x.Inverse(&x)
		e = new(big.Int).Neg(k)
	}
	z.SetOne()
	b := e.Bytes()
	for i := 0; i < len(b); i++ {
		w := b[i]
		for j := 0; j < 8; j++ {
			z.Square(z)
			if (w & (0b10000000 >> j)) != 0 {
				z.Mul(z, &x)
			}
		}
	}
	return z
}

var bigIntPool = sync.Pool{
	New: func() interface{} {
		return new(big.Int)
	},
}
