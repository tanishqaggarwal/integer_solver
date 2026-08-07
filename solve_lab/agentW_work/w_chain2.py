"""W stage 16b: is the 7181 off-pin break MECHANISM or COLLATERAL?"""
import sys, os, re, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import frameB, model
from collections import deque
PVAL = 115792089237316195423570985008687907853269984665640564039457584007908834671663
d = model.get(); A = d['atom_src']; AV = d['atom_vars']
fr = frameB.Frame([642, 28730, 29854, 31864])
W = json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
v0 = [0]*frameB.NV
for k, val in W.items(): v0[int(k[2:]) if k.startswith('x_') else int(k)] = int(val)
fv = {u: v0[u] for u in fr.free if v0[u] != 0}
st = frameB.State(fr, fv); v = fr.forward(fv)
blocks = {b['E']: b for b in json.load(open('w_blocks4.json'))}
b = blocks[7181]
print('block E=7181  gate L=x_%d -> %s   (dead)' % (b['L'], v[b['L']]))
for k in ('i1','i2','i3','i4','i5','i6'):
    x = b[k]; a = fr.definer[x]
    print('  %s = x_%-6d  free=%-5s  def=%s   val%%P==0: %s' %
          (k, x, x in fr.free, (A[a][:70] if a >= 0 else None), v[x] % PVAL == 0))
print()
# where do 33469's equal inputs come from?
d3 = blocks[33469]
for k in ('i1','i2','i3','i4'):
    x = d3[k]
    print('  33469.%s = x_%-6d free=%s ; check atoms it moves: %s' %
          (k, x, x in fr.free, sorted(fr.chk.get(x, []))[:10]))
print()
# do 7181's outputs and 33469's inputs share equations?
eqs = {}
for a in (35759, 35761):
    eqs[a] = set(fr.eq_of[a])
print('equations of the two broken off-pins:', sorted(eqs[35759] | eqs[35761]))
ch = set()
for k in ('i1','i2','i3','i4'):
    ch |= set(fr.chk.get(d3[k], []))
print('check atoms movable by 33469 inputs:', len(ch), sorted(ch)[:16])
e3 = set()
for a in ch: e3 |= set(fr.eq_of[a])
print('equations touched by 33469 inputs:', len(e3))
print('overlap with the off-pin equations:', sorted(e3 & (eqs[35759] | eqs[35761])))
print('overlap with the 7 failing equations:', sorted(e3 & {12231,12270,12350,14584,18673,22044,29125}))
