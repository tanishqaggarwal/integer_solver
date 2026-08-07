"""W stage 11: END-TO-END.  Drive the 39,026 deliverable through the full forward map and
evaluate, EXACTLY OVER Z, every one of the 383 blocks: gate L, A, B, N1, N2, and each of the
3 congruence atoms + 2 off-pins.  Which family is each block actually in?"""
import sys, os, json, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB, model
from collections import Counter
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
d = model.get(); A = d['atom_src']
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st = frameB.State(fr, fv)
print('score from cold:', st.score(), 'failing atoms:', sorted(st.nz()))
v = fr.forward(fv)
blocks = json.load(open('w_blocks4.json'))
ACODE = frameB.ACODE
cls = Counter(); anom = []
detail = []
for b in blocks:
    L = v[b['L']]
    i = [v[b['i%d' % k]] for k in range(1, 7)]
    Q = v[24453]
    Av = i[0]-i[1]; Bv = i[3]-i[2]; Ev = i[0]+i[1]+i[4]+Q
    N1 = Ev*Av*Av - Bv*Bv; N2 = Av*(i[2]+i[5]) - Bv*(i[1]-i[4])
    atoms = [cg['ring']['atom'] for cg in b['congs']]
    vals = [eval(ACODE[a], {'v': v, '__builtins__': {}}) for a in atoms]
    live = (L != 0)
    lawP = (N1 % PVAL == 0 and N2 % PVAL == 0)
    deg = (Av % PVAL == 0 and Bv % PVAL == 0)
    zero_in = all(x == 0 for x in i[:4])
    tag = ('DEAD(all inputs 0)' if zero_in else
           ('gate off' if not live else
            ('DEGENERACY A=B=0 mod P' if deg else ('chord' if lawP else 'LAW VIOLATED'))))
    cls[(tag, tuple(x == 0 for x in vals))] += 1
    if tag == 'LAW VIOLATED' or any(x != 0 for x in vals):
        anom.append((b['E'], tag, L, [x == 0 for x in vals], N1 % PVAL == 0, N2 % PVAL == 0))
    detail.append(dict(E=b['E'], L=int(L), tag=tag, atomsz=[x == 0 for x in vals],
                       N1z=N1 % PVAL == 0, N2z=N2 % PVAL == 0,
                       Az=Av % PVAL == 0, Bz=Bv % PVAL == 0))
print()
for k, n in sorted(cls.items(), key=lambda t: -t[1]): print('  %4d  %s  congruence-atoms-zero=%s' % (n, k[0], k[1]))
print()
print('anomalous blocks:', len(anom))
for t in anom[:12]: print('   ', t)
json.dump(detail, open('w_deliv.json', 'w'))
