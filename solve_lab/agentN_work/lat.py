"""Two-knob lattice optimisation of the residual region in the agent-H frame."""
import ev, fast, json, os, time, itertools, sys
from fast import St
from close2 import close
p=115792089237316195423570985008687907853269984665640564039457584007908834671663
C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
BITS=json.load(open('bits.json')); ALL=set(BITS['A']+BITS['B'])
b=int(sys.argv[1]) if len(sys.argv)>1 else 47
st0=St({})
st=st0.clone().set_free({b:1})
out,ok,fr=close(st,frozen=set(ALL),maxsteps=300)
print('base score',out.score(),'nz',sorted(out.nz()),'fails',sorted(out.fails))
base=out
def shift(st,s,t):
    g=st.clone()
    g.set_free({18956: st.fv.get(18956,0)+p*s, 7497: st.fv.get(7497,0)+8863713*s,
                24468: st.fv.get(24468,0)+p*t, 11436: st.fv.get(11436,0)+t})
    return g
g00=base; g10=shift(base,1,0); g01=shift(base,0,1)
print('shift(1,0) score',g10.score(),'nz',sorted(g10.nz()))
print('shift(0,1) score',g01.score(),'nz',sorted(g01.nz()))
E=sorted(set(g00.fails)|set(g10.fails)|set(g01.fails))
print('region equations',len(E))
rows=[]
for e in E:
    c=g00.eq[e]; a1=g10.eq[e]-c; a2=g01.eq[e]-c
    rows.append((e,a1,a2,c))
    # verify linearity
g11=shift(base,1,1)
lin=all(g11.eq[e]==c+a1+a2 for e,a1,a2,c in rows)
print('linear in (s,t):',lin)
for e,a1,a2,c in rows[:30]:
    print(' eq%d a=%s b=%s c=%s'%(e,a1,a2,c))
json.dump([[e,str(a1),str(a2),str(c)] for e,a1,a2,c in rows],open('lat_rows_%d.json'%b,'w'))
