"""Incremental exact evaluator: change a few free inputs, update only affected cone."""
import ev, pickle, os, json
from collections import defaultdict
HERE=os.path.dirname(os.path.abspath(__file__))
S=pickle.load(open(os.path.join(HERE,'support.pkl'),'rb'))
sup=S['sup']; csup=S['csup']; fidx=S['fidx']; inv={i:v for v,i in fidx.items()}
NV=38748
order=ev.F['order']; definer=ev.F['definer']
pos={v:i for i,v in enumerate(order)}          # topological position of defined vars
desc=defaultdict(list)                          # free input -> defined vars depending on it
for v in order:
    s=sup[v]
    while s:
        b=s&-s; desc[inv[b.bit_length()-1]].append(v); s^=b
for u in desc: desc[u].sort(key=lambda v:pos[v])
chk=defaultdict(list)                            # free input -> check atoms depending on it
for a,s in csup.items():
    t=s
    while t:
        b=t&-t; chk[inv[b.bit_length()-1]].append(a); t^=b
eq_of=defaultdict(list)
for i,(m,sq,tl) in enumerate(ev.eq_terms):
    for c,a in tl: eq_of[a].append(i)
DEF={v:c for v,c in ev.DEFSEQ}
NEQ=len(ev.eq_terms)

class St:
    __slots__=('fv','v','av','eq','fails','ns')
    def __init__(s,fv):
        s.fv=dict(fv)
        s.v=ev.forward(s.fv)
        s.ns={'v':s.v,'__builtins__':{}}
        s.av={a:eval(c,s.ns) for a,c in ev.CHECKCODE.items()}
        s.eq=[0]*NEQ
        for i,(m,sq,tl) in enumerate(ev.eq_terms):
            t=0
            for c,a in tl:
                x=s.av.get(a)
                if x: t+=c*x
            s.eq[i]=m*(t*t if sq else t)
        s.fails=set(i for i in range(NEQ) if s.eq[i])
    def clone(s):
        o=St.__new__(St)
        o.fv=dict(s.fv); o.v=s.v[:]; o.ns={'v':o.v,'__builtins__':{}}
        o.av=dict(s.av); o.eq=s.eq[:]; o.fails=set(s.fails)
        return o
    def set_free(s,changes):
        """changes: dict free_input -> value.  Returns self (mutated)."""
        aff=set(); ch=set()
        for u,val in changes.items():
            if s.v[u]==val and s.fv.get(u,0)==val: continue
            s.fv[u]=val; s.v[u]=val
            aff.update(desc[u]); ch.update(chk[u])
        for v in sorted(aff,key=lambda v:pos[v]):
            s.v[v]=eval(DEF[v],s.ns)
        eqs=set()
        for a in ch:
            nv=eval(ev.CHECKCODE[a],s.ns)
            if nv!=s.av[a]:
                s.av[a]=nv; eqs.update(eq_of[a])
        for i in eqs:
            m,sq,tl=ev.eq_terms[i]
            t=0
            for c,a in tl:
                x=s.av.get(a)
                if x: t+=c*x
            val=m*(t*t if sq else t)
            s.eq[i]=val
            if val: s.fails.add(i)
            else: s.fails.discard(i)
        return s
    def score(s): return NEQ-len(s.fails)
    def nz(s): return [a for a,x in s.av.items() if x]

if __name__=='__main__':
    import time
    t0=time.time(); st=St({}); print('base score',st.score(),'nz',sorted(st.nz()),'%.1fs'%(time.time()-t0))
    C1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
    C2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
    t0=time.time()
    st2=st.clone().set_free({5096:C2,21589:C1,16742:C1,12186:C2,18956:C1,24468:C2,542:1,1438:1})
    print('mux score',st2.score(),'nz',sorted(st2.nz()),'%.3fs'%(time.time()-t0))
    t0=time.time()
    for _ in range(200):
        st3=st2.clone().set_free({13153:12345})
    print('200 clone+update %.2fs'%(time.time()-t0))
