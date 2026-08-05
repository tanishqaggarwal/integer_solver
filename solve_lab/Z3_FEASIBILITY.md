# z3/SMT feasibility for the trapdoor residual

z3 5.0.0 available. Measured (Int + explicit quotient encoding, p=2^256-2^32-977):
- Single modular multiply `x*y ≡ c mod p`: SAT in 0.0s.
- Single modular sqrt `z^2 ≡ a mod p` (a a QR): SAT in 0.0s.
- **Two coupled products `x*y≡A AND (x+y)^2≡B mod p`: UNKNOWN (timeout 60s).**

Conclusion: z3 handles ISOLATED 256-bit modular constraints instantly but chokes
the moment products couple — which is exactly the value×value core / CRT chain
that walls every approach. So SMT is NOT a black-box solver here; its only viable
use is the modular/CRT feasibility layer with a CEGAR loop (gate outputs as vars →
linear-mod-p atoms → forward-check realizability → add one linearized cut), never
the full coupled product web in a single call. Guidance relayed to the SMT agent.
