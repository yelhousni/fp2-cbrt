# fp2-cbrt

Companion code for *"Fast cube roots in Fp2 via the algebraic torus"*.

Implements three methods for computing cube roots in quadratic extension fields Fp2:

- **Torus**: reduces the Fp2 cube root to base-field Fp operations via the algebraic torus T2(Fp) and Lucas V-sequences (Algorithm 1 in the paper).
- **Frobenius**: Frobenius-splitting multi-exponentiation in Fp2 using Strauss-Shamir's trick.
- **Direct**: standard exponentiation in Fp2.

## Target primes

| Name | Prime | log2(p) | p mod 9 | Source |
|------|-------|---------|---------|--------|
| p-p254 | BN254 base field | 254 | 1 | Pairing |
| p-p377 | BLS12-377 base field | 377 | 7 | Pairing |
| p-p381 | BLS12-381 base field | 381 | 1 | Pairing |
| i-p381 | 2^372 * 437 - 1 | 381 | 4 | Isogeny |
| i-p575 | 2^567 * 139 - 1 | 575 | 4 | Isogeny |
| i-p765 | 2^756 * 257 - 1 | 765 | 4 | Isogeny |

## Repository structure

```
fields/          Pre-generated Fp arithmetic packages (from gnark-crypto)
  pp254/         BN254 base field
  pp377/         BLS12-377 base field
  pp381/         BLS12-381 base field
  ip381/         2^372 * 437 - 1
  ip575/         2^567 * 139 - 1
  ip765/         2^756 * 257 - 1

fp2/             Per-prime Fp2 arithmetic + cube root implementations
  pp254/
    e2.go        Fp2 arithmetic (Add, Sub, Mul, Square, Inverse, ...)
    cbrt.go      Torus, Frobenius, and Direct cube root methods
    cbrt_test.go Tests and benchmarks
  pp377/
  pp381/
  ip381/
  ip575/
  ip765/
```

## Prerequisites

- Go >= 1.22
- The `fields/` packages depend on `gnark-crypto`. The `go.mod` uses a local `replace` directive pointing to a local gnark-crypto checkout. Adjust the path in `go.mod` if needed:
  ```
  replace github.com/consensys/gnark-crypto => /path/to/gnark-crypto
  ```

## Build & test

```bash
# Run all tests
go test ./fp2/...

# Run tests for a specific prime
go test ./fp2/pp377/

# Verbose
go test -v ./fp2/...
```

## Benchmarks

```bash
# Benchmark all primes
go test -bench=BenchmarkCbrt -benchmem ./fp2/...

# Benchmark a specific prime
go test -bench=BenchmarkCbrt -benchmem -count=4 ./fp2/pp377/

# All three methods for all primes (recommended for paper results)
go test -bench=BenchmarkCbrt -benchmem -count=4 ./fp2/...
```

### Sample results (Apple M5, ARM64, Go 1.24)

| Prime | Direct (us) | Frobenius (us) | Torus (us) | vs Direct | vs Frobenius |
|-------|------------|----------------|-----------|-----------|--------------|
| p-p254 | 20.6 | 10.5 | 9.7 | 2.13x | 1.08x |
| p-p377 | 69.9 | 34.1 | 25.0 | 2.80x | 1.37x |
| p-p381 | 60.6 | 29.1 | 25.1 | 2.41x | 1.16x |
| i-p381 | 47.8 | 24.4 | 19.5 | 2.45x | 1.25x |
| i-p575 | 180.4 | 90.6 | 68.2 | 2.65x | 1.33x |
| i-p765 | 436.0 | 221.1 | 189.2 | 2.30x | 1.17x |

## License

See [LICENSE](LICENSE).
