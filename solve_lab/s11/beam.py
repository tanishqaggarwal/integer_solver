"""Beam search over mirror repairs.

At each state every broken atom is a check that some FREE variable can zero exactly (the
'mirror' shape).  Applying one and rippling moves the defect somewhere else; the response is
often quadratic, so a linear model cannot see it -- but an exact re-base can.  So: beam search
over (broken atom, free variable) moves, scored by the number of failing EQUATIONS, which is
what the objective actually counts.
"""
import sys, os, json, time, heapq, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L
from ip7 import load_raw
P=L.P; HERE=os.path.dirname(os.path.abspath(__file__)); LAB=os.path.join(HERE,'..')
sys.set_int_max_str_digits(400000)
FREE=set(u for u in range(L.NVARS) if u not in L.definer)

def lin_solve(a,t,v):
    c=0
    for m,cc in L.polys[a].items():
        k=m.count(t)
        if k==0: continue
        if k>1: return None
        term=cc
        for u in m:
            if u!=t: term*=v[u]
        c+=term
    if c==0: return None
    old=v[t]; v[t]=0; rest=L.evalpoly(L.polys[a],v); v[t]=old
    if rest%c: return None
    return -rest//c

def brk(v): return [a for a in range(L.NA) if L.evalpoly(L.polys[a],v)!=0]
def nfail(B): return len(set().union(*[set(L.atom2eq.get(a,{})) for a in B])) if B else 0

def moves(v,B):
    out=[]
    for a in B:
        for t in sorted(L.avars[a]):
            if t not in FREE: continue
            x=lin_solve(a,t,v)
            if x is None or x==v[t]: continue
            out.append((a,t,x))
    return out

def search(v0, width=8, depth=25, tlimit=2400):
    B0=brk(v0); f0=nfail(B0)
    beam=[(f0,tuple(sorted(B0)),v0,())]
    best=(f0,v0,())
    seen={tuple(sorted(B0))}
    t0=time.time()
    for d in range(depth):
        nxt=[]
        for f,key,v,path in beam:
            for (a,t,x) in moves(v,list(key)):
                v2=list(v); L.ripple(v2,{t:x})
                B2=brk(v2); f2=nfail(B2); k2=tuple(sorted(B2))
                if k2 in seen: continue
                seen.add(k2)
                nxt.append((f2,k2,v2,path+((a,t),)))
                if f2<best[0]:
                    best=(f2,v2,path+((a,t),))
                    print(f"   depth{d}: NEW BEST failing={f2} score={L.NEQ-f2} atoms={list(k2)}", flush=True)
        if not nxt: break
        nxt.sort(key=lambda z:(z[0],len(z[1])))
        beam=nxt[:width]
        print(f"  depth{d}: {len(nxt)} states, beam best={beam[0][0]} atoms={list(beam[0][1])[:10]} ({time.time()-t0:.0f}s)", flush=True)
        if time.time()-t0>tlimit: break
    return best

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'data','fix2_round.json')
    v=load_raw(src)
    f,vb,path=search(v)
    B=brk(vb)
    print(f"BEST failing={f} score={L.NEQ-f} broken atoms={B}")
    print("path:",path)
    json.dump({('x_%d'%i):vb[i] for i in range(L.NVARS)}, open(os.path.join(HERE,'data','beam_out.json'),'w'))
