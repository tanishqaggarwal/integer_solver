"""Corrected cut family: inject the forged value AND propagate it up the branch to m."""
import sys, os, json, pickle, time, collections
from math import gcd
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import checker as CK
src=open('/home/user/integer_solver/solve_lab/agentL_work/cansearch.py').read().split("print('loading checker...'")[0]
src='\n'.join(l for l in src.split('\n') if 'import checker as CK' not in l)
exec(src)
print('loading checker...',flush=True); t0=time.time()
CODES,_=CK.load_equations(); print('  %.0fs'%(time.time()-t0),flush=True)
def exact_fail(vv):
    v=[0]*CK.NVARS; n=min(len(vv),CK.NVARS); v[:n]=vv[:n]
    ns={'v':v,'__builtins__':{}}
    return [i for i,c in enumerate(CODES) if eval(c,ns)!=0]
def build2(G,L,csite,vabmode='zero',rounds=60):
    m=LCA(G,L); S={G,L}
    v,isl,valn=assignment(S,ORIENT); v[24468]=T1; v[18956]=T2
    ga=NODE[m]['a']; gb=NODE[m]['b']
    chG,chL=(ga,gb) if (G in sub[ga]) else (gb,ga)
    sG='va' if chG==ga else 'vb'
    pmG=perm[(m,sG)]
    GV=(valn[chG][pmG[0]],valn[chG][pmG[1]])           # G's value, in m's frame
    # inject at csite AND propagate all the way up to m
    y=csite
    while y!=m:
        pn=parent[y]; ps=side_of[y]
        Vp=down(GV,m,pn) if pn!=m else GV              # forged value in pn's own frame
        for i,d in enumerate(OUT[pn]): v[d[ps]]=Vp[i]
        y=pn
    if vabmode=='forge':
        W=down(GV,m,csite)
        for i,d in enumerate(OUT[csite]): v[d['vab']]=W[i]
    elif vabmode=='deliv':
        D=json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json'))
        for d in OUT[csite]:
            v[d['vab']]=int(D.get('x_%d'%d['vab'],0))
    # m now sees two equal inputs -> its output is unconstrained; drive it to TARGET
    req=down(TGT,ROOT,m) if m!=ROOT else TGT
    for i,d in enumerate(OUT[m]): v[d['vab']]=req[i]
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    for rd in range(rounds):
        bad=relift(vv)
        if not bad: break
        r=E.run(vv); fixed=0
        for a in bad:
            i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
            if cur%p: continue
            imm=[q for q in vars_of(E.atoms[a]) if q in SHIFT]
            for w in imm+[q for q in atomvalvars[a] if q in SHIFT and q not in imm]:
                old=vv[w]; vv[w]=old+p; d2=E.run(vv)[i]-cur; vv[w]=old
                if d2==0: continue
                g=gcd(d2,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d2//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; fixed+=1; break
        if fixed==0: break
    relift(vv)
    return vv
if __name__=='__main__':
    for mode in ('deliv','forge','zero'):
        vv=build2(24601,2081,27994,vabmode=mode)
        r=E.run(vv); nz=[E.res[i] for i,x in enumerate(r) if x]
        f=exact_fail(vv)
        print('site=27994 vab=%-6s -> %d nonzero atoms, EXACT failing %d %s'%(mode,len(nz),len(f),sorted(f)[:10]),flush=True)
        for a in nz: print('        ',a[:120])
        if len(f)<7:
            json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('cut_%s.json'%mode,'w'))
            print('        *** WROTE cut_%s.json  (%d failing)'%(mode,len(f)),flush=True)
