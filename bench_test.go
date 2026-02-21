package main

import (
	"testing"

	e2ip381 "github.com/yelhousni/fp2-cbrt/fp2/ip381"
	e2ip575 "github.com/yelhousni/fp2-cbrt/fp2/ip575"
	e2pp254 "github.com/yelhousni/fp2-cbrt/fp2/pp254"
	e2pp377 "github.com/yelhousni/fp2-cbrt/fp2/pp377"
	e2pp381 "github.com/yelhousni/fp2-cbrt/fp2/pp381"
)

func BenchmarkCbrt(b *testing.B) {
	b.Run("pp254/Torus", func(b *testing.B) {
		var a, x e2pp254.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.Cbrt(&x)
		}
	})
	b.Run("pp254/Direct", func(b *testing.B) {
		var a, x e2pp254.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.CbrtDirect(&x)
		}
	})

	b.Run("pp381/Torus", func(b *testing.B) {
		var a, x e2pp381.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.Cbrt(&x)
		}
	})
	b.Run("pp381/Direct", func(b *testing.B) {
		var a, x e2pp381.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.CbrtDirect(&x)
		}
	})

	b.Run("pp377/Torus", func(b *testing.B) {
		var a, x e2pp377.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.Cbrt(&x)
		}
	})
	b.Run("pp377/Direct", func(b *testing.B) {
		var a, x e2pp377.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.CbrtDirect(&x)
		}
	})

	b.Run("ip381/Torus", func(b *testing.B) {
		var a, x e2ip381.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.Cbrt(&x)
		}
	})
	b.Run("ip381/Direct", func(b *testing.B) {
		var a, x e2ip381.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.CbrtDirect(&x)
		}
	})

	b.Run("ip575/Torus", func(b *testing.B) {
		var a, x e2ip575.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.Cbrt(&x)
		}
	})
	b.Run("ip575/Direct", func(b *testing.B) {
		var a, x e2ip575.E2
		a.SetRandom()
		x.Square(&a)
		x.Mul(&x, &a)
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			a.CbrtDirect(&x)
		}
	})
}
