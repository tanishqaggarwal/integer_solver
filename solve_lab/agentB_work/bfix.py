"""Agent B: drive the four residual divisibility conditions using the free inputs
x4432 (-> x28730), x14865 (-> x8731), x31861 (-> x9118), x7068 (-> x7068-x2099).
Iterate: forward-evaluate, measure residues, correct, repeat."""
import os, sys, json, pickle, collections, time
os.environ.setdefault('ORIENT', 'orient3.pkl')
sys.path.insert(0, '/home/user/integer_solver/solve_lab/agentB_work')
import beval as E

P = 115792089237316195423570985008687907853269984665640564039457584007908834671663
Q = 7376877 * P
NV = 38748
KNOBS = {'x4432': 4432, 'x14865': 14865, 'x31861': 31861, 'x7068': 7068}
HANDLES = [642, 31864, 29854, 9413, 10903, 17325, 1329]

def residues(val):
    return {
        'A p|x28730': val[28730] % P,
        'B p|x8731': val[8731] % P,
        'E p|x9118': val[9118] % P,
        'D 7376877p|(x7068-x2099)': (val[7068] - val[2099]) % Q,
    }

def set_handles(val):
    """Set the residual quotient handles ONLY when the congruence actually holds."""
    if val[8731] % P == 0:
        val[31864] = -val[7075]*val[8731]
        val[10903] = val[31864]//P
    if val[9118] % P == 0:
        val[29854] = 5113045*val[7075]*val[9118]
        val[1329] = val[29854]//P
    if val[28730] % P == 0:
        val[9413] = val[28730]//P
    d = val[7068] - val[2099]
    if d % (7376877*P) == 0:
        val[642] = d//7376877
        val[17325] = val[642]//P
    return val

def run(base, deltas, iters=6, verbose=True):
    v0 = E.load(base)
    freeval = {v: v0[v] for v in E.free}
    for k, d in deltas.items(): freeval[k] = freeval.get(k, v0[k]) + d
    val, nd, nf2 = E.forward(freeval, default=v0)
    for it in range(iters):
        r = residues(val)
        if verbose: print("  iter %d residues: %s" % (it, {k: (v.bit_length() if v else 0) for k, v in r.items()}))
        if all(x == 0 for x in r.values()): break
        # corrections
        freeval[4432] = freeval.get(4432, v0[4432]) - (val[28730] % P)
        freeval[14865] = freeval.get(14865, v0[14865]) + (val[8731] % P)
        freeval[31861] = freeval.get(31861, v0[31861]) + (val[9118] % P)
        freeval[7068] = freeval.get(7068, v0[7068]) - ((val[7068]-val[2099]) % Q)
        val, nd, nf2 = E.forward(freeval, default=v0)
    val = set_handles(val)
    return val, v0

if __name__ == '__main__':
    base = '/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'
    val, v0 = run(base, {})
    print("final residues:", {k: (v.bit_length() if v else 0) for k, v in residues(val).items()})
    ok, fail, fv = E.score(val)
    print("SCORE:", ok, fail[:15])
    out = sys.argv[1] if len(sys.argv) > 1 else 'out/bfix.json'
    json.dump({('x_%d'%i): val[i] for i in range(NV)}, open(out, 'w'))
    print("wrote", out)
