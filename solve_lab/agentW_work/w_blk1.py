"""W: locate the 383 law blocks from EQUATIONS.txt via my own model, no other agent's code."""
import sys, os, json, re, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.set_int_max_str_digits(20_000_000)
import model
d = model.get(); A_SRC=d['atom_src']; A_VAR=d['atom_vars']; EQ=d['eq_terms']
print('atoms',len(A_SRC),'eqs',len(EQ))
# constants: find atoms of the form x_v - CONST with CONST huge
CONST={}
pat=re.compile(r'^\(*x_(\d+)\)*\s*-\s*\(*(-?\d+)\)*$')
big=0
for i,s in enumerate(A_SRC):
    t=s.replace(' ','')
    m=re.match(r'^\(*x_(\d+)\)*-\(+(-?\d+)\)+$', t)
    if m:
        v=int(m.group(1)); c=int(m.group(2))
        CONST.setdefault(v,set()).add(c)
        if abs(c)>2**200: big+=1
print('pinned vars',len(CONST),'big pins',big)
bigv={v:list(cs)[0] for v,cs in CONST.items() if len(cs)==1 and abs(list(cs)[0])>2**200}
print('vars pinned to >200-bit const:',len(bigv))
from collections import Counter
cnt=Counter(bigv.values())
print('distinct big consts:',len(cnt))
# how many atoms reference each big var
ref=Counter()
for i,vs in enumerate(A_VAR):
    for v in vs:
        if v in bigv: ref[v]+=1
print('top referenced big-const vars:', [(v,ref[v],bigv[v].bit_length()) for v,_ in ref.most_common(6)])
json.dump({str(k):str(v) for k,v in bigv.items()}, open('w_bigconst.json','w'))
