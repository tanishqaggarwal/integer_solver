"""W stage 18: THE OUT-OF-K MEASUREMENT.
At block 7181 the deliverable injects (i5,i6 nonzero mod P at a DEAD gate) and pays 5 of the
site's 9 equations.  Do the same at the four equation-disjoint minimum-incidence blocks and
price it EXACTLY through frameB.State -- no linear model."""
import sys, os, json, itertools, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
A = frameB.atom_src
B = {b['E']: b for b in json.load(open('w_blocks4.json'))}
PRICE = {r['E']: r for r in json.load(open('w_price.json'))}
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)

fr = frameB.Frame([642, 28730, 29854, 31864])
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st0 = frameB.State(fr, fv)
print('frame B base:', st0.score(), sorted(st0.fails))
print()
print('--- control: what does ZEROING block 7181 s outputs cost? ---')
for comb in ([9118], [8731], [9118, 8731]):
    s = st0.clone().set_free({u: 0 for u in comb})
    print('   zero %-16s -> score %d  fails %s' % (comb, s.score(), sorted(s.fails)))
print()
print('--- the five minimum-incidence blocks: are their slots free inputs? gate value? ---')
vdef = fr.forward(fv)
for E in (7181, 3227, 4429, 30886, 31606):
    b = B[E]
    sl = [b['i%d' % k] for k in range(1, 7)]
    print('  E=%-6d gate x_%-6d = %-3s  slots free: %s' %
          (E, b['L'], vdef[b['L']], [(x, x in fr.free, vdef[x] % PVAL == 0) for x in sl]))
    print('        9 equations: %s' % PRICE[E]['neq'], sorted({e for a in PRICE[E]['offpins'] for e in fr.eq_of[a]}))
