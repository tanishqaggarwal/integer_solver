"""Cut at node x: put the target-inverted value on x's parent slot wires. Only those 2 atoms break."""
import sys, pickle, json, collections
src=open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/mkassign2.py').read().split('#MAINSTART')[0]
exec(src)
TGT=tuple(pickle.load(open('target.pkl','rb')))
NODE=M['NODE']; OUT=M['OUT']; tree=M['tree']; sub=M['sub']; ROOT=M['ROOT']; liveset=set(M['live'])
parent={}; side_of={}
for n in NODE:
    for s,ch in (('va',NODE[n]['a']),('vb',NODE[n]['b'])): parent[ch]=n; side_of[ch]=s
def cutbuild(x):
    n=parent[x]; side=side_of[x]
    cand=[l for l in sub[x] if l in liveset]
    if not cand: return None
    S={cand[0]}
    v,isl,valn=assignment(S,ORIENT)
    v[24468]=T1; v[18956]=T2
    chain=[]; y=n
    while y!=ROOT: chain.append(y); y=parent[y]
    chain.append(ROOT); chain.reverse()
    req=TGT
    for k,m in enumerate(chain):
        pc = chain[k+1] if k+1<len(chain) else x
        sk = 'va' if NODE[m]['a']==pc else 'vb'
        other = NODE[m]['b'] if sk=='va' else NODE[m]['a']
        if isl.get(other): return None          # only the all-siblings-dead case here
        R = req
        if m==n:
            for i,d in enumerate(OUT[n]): v[d[side]] = R[i]
            break
        pm=perm[(m,sk)]; req=(R[pm[0]],R[pm[1]])
    vv=[0]*NV
    for k2,xx in v.items(): vv[k2]=xx
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
    r=E.run(vv)
    nz=[i for i,q in enumerate(r) if q]
    return vv,nz,E.score(r)
if __name__=='__main__':
    PR=pickle.load(open('price.pkl','rb'))
    best=None
    for cost,n,side,ws in PR['rows'][:40]:
        x = NODE[n]['a'] if side=='va' else NODE[n]['b']
        out=cutbuild(x)
        if out is None: continue
        vv,nz,bad=out
        print('node x%-6d %s wires %s  incidence-cost %-3d  atoms nonzero %d  F-model failing %d'%(
            n,side,ws,cost,len(nz),len(bad)),flush=True)
        if best is None or len(bad)<best[0]:
            best=(len(bad),n,side,vv)
            json.dump({'x_%d'%i:vv[i] for i in range(NV) if vv[i]},open('cut_best.json','w'))
    print('best F-model failing:',best[0] if best else None,'node',best[1:3] if best else None)
