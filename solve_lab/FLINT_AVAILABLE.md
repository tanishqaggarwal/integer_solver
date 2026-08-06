# python-flint 0.9.0 installed (pip)

Fast exact 256-bit GF(p) algebra now available:
- `flint.fmpz_mod_mat` — matrices mod p (rank, nullspace, solve) — MUCH faster than pure-Python Gaussian.
- `flint.fmpz_mod_poly` — univariate mod p: `.factor()`, `.roots()` — fast Tonelli-Shanks/Cantor-Zassenhaus.
- `flint.fmpz_mod_mpoly` — multivariate over Z/pZ (256-bit modulus). Context: `fmpz_mod_mpoly_ctx.get(names, p, flint.Ordering.lex)`.
- NO Gröbner over GF(p) (buchberger_naive is only on integer fmpz_mpoly_vec, not fmpz_mod). XL infeasible: residual support = 256 bits.

Use for: fast linear-algebra steps (null space, XL degree-2), modular root-finding in core/cube solves.
