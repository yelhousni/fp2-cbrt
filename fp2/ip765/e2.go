// Fp2 = Fp[u]/(u² + 1) arithmetic for isogeny prime 2^756·257-1
package ip765

import (
	"math/big"
	"sync"

	fp "github.com/yelhousni/fp2-cbrt/fp/ip765"
)

type E2 struct {
	A0, A1 fp.Element
}

func (z *E2) SetZero() *E2     { z.A0.SetZero(); z.A1.SetZero(); return z }
func (z *E2) SetOne() *E2      { z.A0.SetOne(); z.A1.SetZero(); return z }
func (z *E2) Set(x *E2) *E2    { z.A0.Set(&x.A0); z.A1.Set(&x.A1); return z }
func (z *E2) Equal(x *E2) bool { return z.A0.Equal(&x.A0) && z.A1.Equal(&x.A1) }
func (z *E2) IsZero() bool     { return z.A0.IsZero() && z.A1.IsZero() }

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

// Mul: β = -1
func (z *E2) Mul(x, y *E2) *E2 {
	var a, b, c fp.Element
	a.Add(&x.A0, &x.A1)
	b.Add(&y.A0, &y.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &y.A0)
	c.Mul(&x.A1, &y.A1)
	z.A1.Sub(&a, &b).Sub(&z.A1, &c)
	z.A0.Sub(&b, &c)
	return z
}

func (z *E2) Square(x *E2) *E2 {
	var a, b fp.Element
	a.Add(&x.A0, &x.A1)
	b.Sub(&x.A0, &x.A1)
	a.Mul(&a, &b)
	b.Mul(&x.A0, &x.A1).Double(&b)
	z.A0.Set(&a)
	z.A1.Set(&b)
	return z
}

func (z *E2) Inverse(x *E2) *E2 {
	var t0, t1 fp.Element
	t0.Square(&x.A0)
	t1.Square(&x.A1)
	t0.Add(&t0, &t1)
	t1.Inverse(&t0)
	z.A0.Mul(&x.A0, &t1)
	z.A1.Mul(&x.A1, &t1).Neg(&z.A1)
	return z
}

func (z *E2) norm(x *fp.Element) {
	var tmp fp.Element
	x.Square(&z.A0)
	tmp.Square(&z.A1)
	x.Add(x, &tmp)
}

func (z *E2) MulByElement(x *E2, y *fp.Element) *E2 {
	var yCopy fp.Element
	yCopy.Set(y)
	z.A0.Mul(&x.A0, &yCopy)
	z.A1.Mul(&x.A1, &yCopy)
	return z
}

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
	New: func() interface{} { return new(big.Int) },
}
