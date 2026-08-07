"""Cancellation search over the deliverable's ACTUAL cut family, priced with the exact checker."""
import sys, os, json, pickle, time, collections
from math import gcd
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import checker as CK
src=open('/home/user/integer_solver/solve_lab/agentL_work/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
TGT=tuple(pickle.load(open('target.pkl','rb')))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; sub=M['sub']; ROOT=M['ROOT']; liveset=set(M['live'])
parent={}; side_of={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=s
def anc(x):
    r=[]
    while x!=ROOT: x=parent[x]; r.append(x)
    return r
def LCA(a,b):
    sa=set(anc(a)+[a])
    x=b
    while x not in sa: x=parent[x]
    return x
def down(val, frm, to):
    """transport a value from ancestor frame `frm` down to node `to`'s own frame"""
    ch=[]; x=to
    while x!=frm: ch.append(x); x=parent[x]
    for y in reversed(ch):
        pm=perm[(parent[y],side_of[y])]
        val=(val[pm[0]],val[pm[1]])
    return val
print('loading checker...',flush=True); t0=time.time()
CODES,_=CK.load_equations(); print('  %.0fs'%(time.time()-t0),flush=True)
def exact_fail(vv):
    v=[0]*CK.NVARS; n=min(len(vv),CK.NVARS); v[:n]=vv[:n]
    ns={'v':v,'__builtins__':{}}
    return sum(1 for c in CODES if eval(c,ns)!=0)
def build(G,L,csite,set_vab=True):
    """ON={G,L}; overwrite the L-branch at csite so both inputs of m=LCA(G,L) are equal;
       then m's output is free -> set it so everything above lands on TARGET."""
    m=LCA(G,L)
    if csite is None or m in set(anc(csite))^set() and False: pass
    S={G,L}
    v,isl,valn=assignment(S,ORIENT); v[24468]=T1; v[18956]=T2
    ga=NODE[m]['a']; gb=NODE[m]['b']
    chG,chL=(ga,gb) if (G in sub[ga]) else (gb,ga)
    sG='va' if chG==ga else 'vb'; sL='vb' if sG=='va' else 'va'
    pmG=perm[(m,sG)]
    GV=(valn[chG][pmG[0]],valn[chG][pmG[1]])          # G's value in m's frame
    # 1. overwrite the L branch at csite (csite is on L's path, strictly under m)
    W=down(GV,m,csite)                                 # value in csite's own frame
    pn=parent[csite]; ps=side_of[csite]
    pmc=perm[(pn,ps)]
    for i,d in enumerate(OUT[pn]): v[d[ps]]=W[pmc[i]]   # parent's slot wires (breaks 2 link atoms)
    if set_vab:
        for i,d in enumerate(OUT[csite]): v[d['vab']]=W[i]   # breaks 2 guard atoms
    # 2. m now sees equal inputs -> its output is unconstrained; drive it to TARGET
    req=down(TGT,ROOT,m) if m!=ROOT else TGT
    for i,d in enumerate(OUT[m]): v[d['vab']]=req[i]
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    for rd in range(60):
        bad=relift(vv)
        if not bad: break
        r=E.run(vv); fixed=0
        for a in bad:
            i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
            if cur%p: continue
            imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
                if d==0: continue
                g=gcd(d,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; fixed+=1; break
        if fixed==0: break
    relift(vv)
    return vv
if __name__=='__main__':
    G,L=24601,2081
    # the deliverable's own site: csite = x27994 (the node holding leaf 2081)
    for sv in (True,False):
        vv=build(G,L,27994,set_vab=sv)
        r=E.run(vv); nz=sum(1 for x in r if x)
        print('VALIDATION G=24601 L=2081 csite=27994 set_vab=%s -> nonzero atoms %d, EXACT failing %d  (deliverable=7)'%(
            sv,nz,exact_fail(vv)),flush=True)
