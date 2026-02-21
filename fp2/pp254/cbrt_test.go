package pp254

import (
	"testing"
)

func TestCbrtTorus(t *testing.T) {
	for i := 0; i < 50; i++ {
		// Generate a random cube: x = a³
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z E2
		res := z.Cbrt(&x)
		if res == nil {
			t.Fatalf("Cbrt returned nil for a cube (iter %d)", i)
		}

		// Verify: z³ = x
		var check E2
		check.Square(&z)
		check.Mul(&check, &z)
		if !check.Equal(&x) {
			t.Fatalf("z³ != x (iter %d)", i)
		}
	}
}

func TestCbrtTorusBaseField(t *testing.T) {
	for i := 0; i < 20; i++ {
		var a E2
		a.A0.SetRandom()
		a.A1.SetZero()
		var x E2
		x.Square(&a)
		x.Mul(&x, &a)

		var z E2
		res := z.Cbrt(&x)
		if res == nil {
			t.Fatalf("Cbrt returned nil for base field cube (iter %d)", i)
		}
		var check E2
		check.Square(&z)
		check.Mul(&check, &z)
		if !check.Equal(&x) {
			t.Fatalf("z³ != x for base field (iter %d)", i)
		}
	}
}

func TestCbrtDirectVsTorus(t *testing.T) {
	for i := 0; i < 10; i++ {
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z1, z2 E2
		r1 := z1.Cbrt(&x)
		r2 := z2.cbrtDirect(&x)
		if r1 == nil || r2 == nil {
			t.Fatalf("one method returned nil (iter %d)", i)
		}
		// Both should be valid cube roots (may differ by cube root of unity)
		var c1, c2 E2
		c1.Square(&z1).Mul(&c1, &z1)
		c2.Square(&z2).Mul(&c2, &z2)
		if !c1.Equal(&x) || !c2.Equal(&x) {
			t.Fatalf("cube root verification failed (iter %d)", i)
		}
	}
}

func TestFp2Arithmetic(t *testing.T) {
	for i := 0; i < 50; i++ {
		var a, b, c E2
		a.SetRandom()
		b.SetRandom()

		// a + b - b = a
		c.Add(&a, &b).Sub(&c, &b)
		if !c.Equal(&a) {
			t.Fatal("add/sub failed")
		}

		// a * b * b⁻¹ = a
		var d E2
		d.Inverse(&b)
		c.Mul(&a, &b).Mul(&c, &d)
		if !c.Equal(&a) {
			t.Fatal("mul/inv failed")
		}

		// square vs mul
		c.Square(&a)
		d.Mul(&a, &a)
		if !c.Equal(&d) {
			t.Fatal("square != mul(a,a)")
		}
	}
}

func BenchmarkCbrtTorus(b *testing.B) {
	var a, x E2
	a.SetRandom()
	x.Square(&a)
	x.Mul(&x, &a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.Cbrt(&x)
	}
}

func BenchmarkCbrtDirect(b *testing.B) {
	var a, x E2
	a.SetRandom()
	x.Square(&a)
	x.Mul(&x, &a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.cbrtDirect(&x)
	}
}

func TestCbrtFrobeniusVsTorus(t *testing.T) {
	for i := 0; i < 10; i++ {
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z1, z2 E2
		r1 := z1.Cbrt(&x)
		r2 := z2.cbrtFrobenius(&x)
		if r1 == nil || r2 == nil {
			t.Fatalf("one method returned nil (iter %d)", i)
		}
		var c1, c2 E2
		c1.Square(&z1).Mul(&c1, &z1)
		c2.Square(&z2).Mul(&c2, &z2)
		if !c1.Equal(&x) || !c2.Equal(&x) {
			t.Fatalf("cube root verification failed (iter %d)", i)
		}
	}
}

func BenchmarkCbrtFrobenius(b *testing.B) {
	var a, x E2
	a.SetRandom()
	x.Square(&a)
	x.Mul(&x, &a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.cbrtFrobenius(&x)
	}
}


func TestCbrtTorusPracVsTorus(t *testing.T) {
	for i := 0; i < 10; i++ {
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z1, z2 E2
		r1 := z1.Cbrt(&x)
		r2 := z2.cbrtTorusPrac(&x)
		if r1 == nil || r2 == nil {
			t.Fatalf("one method returned nil (iter %d)", i)
		}
		var c1, c2 E2
		c1.Square(&z1).Mul(&c1, &z1)
		c2.Square(&z2).Mul(&c2, &z2)
		if !c1.Equal(&x) || !c2.Equal(&x) {
			t.Fatalf("cube root verification failed (iter %d)", i)
		}
	}
}

func TestCbrtFrobenius1bitVsTorus(t *testing.T) {
	for i := 0; i < 10; i++ {
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z1, z2 E2
		r1 := z1.Cbrt(&x)
		r2 := z2.cbrtFrobenius1bit(&x)
		if r1 == nil || r2 == nil {
			t.Fatalf("one method returned nil (iter %d)", i)
		}
		var c1, c2 E2
		c1.Square(&z1).Mul(&c1, &z1)
		c2.Square(&z2).Mul(&c2, &z2)
		if !c1.Equal(&x) || !c2.Equal(&x) {
			t.Fatalf("cube root verification failed (iter %d)", i)
		}
	}
}

func BenchmarkCbrtTorusPrac(b *testing.B) {
	var a, x E2
	a.SetRandom()
	x.Square(&a)
	x.Mul(&x, &a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.cbrtTorusPrac(&x)
	}
}

func BenchmarkCbrtFrobenius1bit(b *testing.B) {
	var a, x E2
	a.SetRandom()
	x.Square(&a)
	x.Mul(&x, &a)
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		a.cbrtFrobenius1bit(&x)
	}
}


func TestCbrtAllMethods(t *testing.T) {
	for i := 0; i < 10; i++ {
		var a, x E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)

		var z [5]E2
		methods := []struct {
			name string
			fn   func(*E2, *E2) *E2
		}{
			{"Cbrt", (*E2).Cbrt},
			{"cbrtDirect", (*E2).cbrtDirect},
			{"cbrtFrobenius", (*E2).cbrtFrobenius},
			{"cbrtFrobenius1bit", (*E2).cbrtFrobenius1bit},
			{"cbrtTorusPrac", (*E2).cbrtTorusPrac},
		}

		for j, m := range methods {
			if m.fn(&z[j], &x) == nil {
				t.Fatalf("%s returned nil (iter %d)", m.name, i)
			}
			var check E2
			check.Square(&z[j]).Mul(&check, &z[j])
			if !check.Equal(&x) {
				t.Fatalf("%s: z³ != x (iter %d)", m.name, i)
			}
		}
	}
}
