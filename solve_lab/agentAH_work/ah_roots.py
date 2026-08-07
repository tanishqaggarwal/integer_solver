"""Exact drop-in replacement for closeS3's `rootset_pp`.

`rootset_pp(coeffs,q,e)` in the fleet's single-wire pass is
    {t for t in range(q**e) if peval(coeffs,t,q**e)==0}
i.e. brute force over the whole residue ring.  927 atoms carry a handle cofactor c>1 and the
largest prime power dividing such a c reaches 16,595,977, so one call can cost >10^7 big-integer
Newton evaluations.  That is the entire reason high-|S| runs appear to "stall".

This returns the SAME SET using agent T's own `t_poly.roots_pp` (Hensel lifting over a proper
prime-field root-find), which `t_close2wj` already uses for exactly this purpose in its two-wire
pass.  Semantics are unchanged; only the enumeration is.  `ah_roots_selftest.py` checks equality
against brute force on the actual instance's coefficient shapes.
"""
import sys
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentT_work')
import t_poly as TP

BRUTE = 3000
STATS = {'fast': 0, 'brute': 0, 'fallback': 0, 'all': 0}

def make(peval):
    def rootset_pp_fast(coeffs, q, e):
        m = q**e
        if m <= BRUTE:
            STATS['brute'] += 1
            return {t for t in range(m) if peval(coeffs, t, m) == 0}
        mo = TP.newton_to_mono(list(coeffs), m)
        if mo is None:                      # k! not invertible mod m -- exact fallback
            STATS['fallback'] += 1
            return {t for t in range(m) if peval(coeffs, t, m) == 0}
        rs = TP.roots_pp(mo, m, q, e, brute=BRUTE)
        if rs == 'ALL':                     # identically zero: same set the brute loop returns
            STATS['all'] += 1
            return set(range(m))
        STATS['fast'] += 1
        return set(rs)
    return rootset_pp_fast
