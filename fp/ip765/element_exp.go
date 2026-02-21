package ip765

import "math/big"

// cbrtHelperQMinus4Div9 is (q-4)/9 as a big.Int, precomputed.
var cbrtHelperQMinus4Div9 *big.Int

func init() {
	cbrtHelperQMinus4Div9, _ = new(big.Int).SetString("10823490389574712357131286038814647262939987695285129911641069167011043154307694764188252795924537348992929329133102371256133248462931513091394906404231324673815780404659245489540423676413867196183067241557756414631604319700269283", 10)
}

// ExpByCbrtHelperQMinus4Div9 sets z to x^((q-4)/9) and returns z.
func (z *Element) ExpByCbrtHelperQMinus4Div9(x Element) *Element {
	return z.Exp(x, cbrtHelperQMinus4Div9)
}
