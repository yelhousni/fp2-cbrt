package pp377

import fp "github.com/yelhousni/fp2-cbrt/fp/pp377"

// pracOps encodes the PRAC differential addition chain for the Lucas V-sequence
// exponent e = 3^{-1} mod (p+1), 376 bits.
// Cost: 689 field ops (vs 750 for binary ladder, 8.1% saving).
// Op codes: 0=SWAP, 1..9=Montgomery PRAC cases, 10=FINAL.
var pracOps = [719]byte{
	3, 10, 3, 10, 3, 10, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0,
	3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 3,
	1, 3, 0, 3, 0, 4, 4, 3, 0, 3, 0, 7, 3, 0, 3, 3, 0, 3, 0, 1,
	3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 0, 3, 1, 3, 0, 3, 0, 3, 0, 3,
	0, 4, 3, 3, 0, 3, 1, 3, 0, 2, 0, 9, 9, 9, 9, 9, 4, 4, 4, 5,
	5, 4, 4, 3, 3, 0, 3, 0, 7, 8, 3, 0, 3, 0, 3, 0, 7, 0, 1, 3,
	0, 3, 3, 3, 0, 3, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 0, 4, 3, 0,
	2, 0, 9, 9, 9, 4, 4, 5, 4, 5, 5, 5, 3, 3, 0, 6, 0, 3, 3, 0,
	3, 0, 5, 5, 3, 3, 0, 2, 0, 9, 4, 4, 4, 3, 3, 3, 0, 3, 0, 3,
	3, 0, 3, 3, 3, 0, 3, 3, 3, 0, 3, 3, 3, 0, 3, 0, 4, 3, 3, 3,
	0, 1, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 4, 3, 0, 3, 0, 3, 3,
	3, 0, 3, 1, 1, 3, 0, 3, 0, 3, 0, 5, 4, 3, 3, 0, 3, 0, 5, 4,
	3, 1, 3, 0, 3, 3, 3, 0, 3, 3, 0, 3, 3, 3, 0, 7, 0, 3, 0, 4,
	5, 5, 3, 3, 0, 4, 3, 0, 3, 0, 3, 0, 4, 3, 3, 3, 0, 3, 0, 5,
	4, 5, 3, 3, 0, 6, 0, 3, 0, 3, 0, 1, 3, 0, 4, 4, 3, 0, 3, 0,
	3, 0, 5, 3, 1, 3, 0, 3, 0, 3, 0, 3, 0, 3, 3, 0, 5, 5, 3, 3,
	0, 3, 3, 0, 3, 0, 3, 0, 3, 3, 3, 0, 3, 0, 3, 0, 5, 3, 3, 0,
	7, 8, 8, 8, 3, 0, 3, 0, 3, 3, 0, 3, 0, 7, 8, 3, 0, 5, 4, 3,
	0, 3, 0, 3, 0, 3, 0, 3, 0, 5, 3, 3, 3, 0, 3, 0, 3, 3, 3, 0,
	3, 0, 3, 0, 7, 6, 3, 0, 5, 3, 3, 3, 0, 6, 8, 3, 0, 5, 3, 3,
	0, 3, 3, 0, 3, 0, 3, 3, 0, 3, 3, 0, 3, 0, 3, 0, 3, 3, 0, 3,
	0, 5, 3, 3, 3, 0, 3, 0, 5, 5, 5, 4, 5, 5, 3, 3, 2, 0, 4, 5,
	3, 3, 0, 3, 0, 3, 2, 0, 4, 4, 5, 4, 3, 3, 0, 3, 0, 5, 5, 3,
	2, 0, 4, 4, 4, 3, 0, 3, 0, 1, 3, 0, 1, 3, 0, 3, 3, 3, 0, 3,
	3, 0, 3, 0, 3, 3, 3, 0, 3, 3, 0, 6, 0, 3, 0, 3, 0, 3, 3, 0,
	5, 4, 3, 0, 3, 0, 6, 3, 0, 3, 0, 3, 0, 1, 3, 0, 4, 3, 0, 3,
	0, 3, 3, 3, 0, 3, 3, 0, 3, 0, 7, 7, 0, 3, 0, 4, 3, 3, 0, 3,
	3, 0, 5, 5, 3, 3, 0, 3, 0, 1, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3,
	0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 3, 0, 3, 0, 3, 0, 5, 3, 3, 0,
	4, 3, 0, 2, 0, 4, 5, 5, 3, 3, 3, 0, 3, 0, 3, 3, 0, 5, 3, 3,
	0, 3, 0, 7, 3, 3, 3, 0, 1, 3, 0, 3, 0, 3, 0, 3, 0, 5, 3, 3,
	0, 3, 0, 3, 3, 3, 0, 3, 0, 3, 3, 3, 0, 3, 0, 3, 0, 3, 0, 5,
	4, 3, 3, 1, 3, 0, 1, 3, 0, 3, 0, 3, 3, 0, 5, 3, 3, 0, 7, 7,
	3, 3, 0, 3, 3, 0, 6, 0, 1, 3, 0, 3, 0, 3, 0, 3, 0, 3, 10,
}

// lucasVPrac computes V_e(alpha) using a PRAC differential addition chain
// (Montgomery 1992). Uses 3 registers: A=V_d, B=V_e, C=V_{|d-e|}.
func lucasVPrac(alpha *fp.Element) fp.Element {
	var A, B, C fp.Element
	var T, T2, T3 fp.Element
	var two fp.Element
	two.SetUint64(2)

	A.Set(alpha)
	B.Set(alpha)
	C.Set(&two)

	for _, op := range pracOps {
		switch op {
		case 0: // SWAP
			A, B = B, A
		case 1: // CASE1: 3 Mul
			T.Mul(&A, &B).Sub(&T, &C)
			T2.Mul(&T, &A).Sub(&T2, &B)
			B.Mul(&T, &B).Sub(&B, &A)
			A.Set(&T2)
		case 2, 4: // CASE2, CASE4: 1 Mul + 1 Sq
			T.Mul(&A, &B).Sub(&T, &C)
			A.Square(&A).Sub(&A, &two)
			B.Set(&T)
		case 3: // CASE3: 1 Mul
			T.Mul(&A, &B).Sub(&T, &C)
			C.Set(&B)
			B.Set(&T)
		case 5: // CASE5: 1 Mul + 1 Sq
			T.Mul(&A, &C).Sub(&T, &B)
			A.Square(&A).Sub(&A, &two)
			C.Set(&T)
		case 6: // CASE6: 3 Mul + 1 Sq
			T.Mul(&A, &B).Sub(&T, &C)
			T2.Square(&A).Sub(&T2, &two)
			T3.Mul(&T2, &A).Sub(&T3, &A)
			A.Set(&T3)
			T3.Mul(&T2, &T).Sub(&T3, &C)
			C.Set(&B)
			B.Set(&T3)
		case 7: // CASE7: 3 Mul + 1 Sq
			T.Mul(&A, &B).Sub(&T, &C)
			T2.Square(&A).Sub(&T2, &two)
			T3.Mul(&T2, &A).Sub(&T3, &A)
			T2.Mul(&T, &A).Sub(&T2, &B)
			A.Set(&T3)
			B.Set(&T2)
		case 8: // CASE8: 3 Mul + 1 Sq
			T.Mul(&A, &B).Sub(&T, &C)
			T2.Square(&A).Sub(&T2, &two)
			T3.Mul(&T2, &A).Sub(&T3, &A)
			C.Mul(&A, &C).Sub(&C, &B)
			A.Set(&T3)
			B.Set(&T)
		case 9: // CASE9: 1 Mul + 1 Sq
			T.Mul(&C, &B).Sub(&T, &A)
			B.Square(&B).Sub(&B, &two)
			C.Set(&T)
		case 10: // FINAL: 1 Mul
			A.Mul(&A, &B).Sub(&A, &C)
			B.Set(&A)
			C.Set(&two)
		}
	}

	return A
}
