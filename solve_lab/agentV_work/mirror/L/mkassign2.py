import sys, pickle, json, time, collections
from math import gcd
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/lift.py').read().split("if __name__")[0]
exec(src)
SL=pickle.load(open('slopes.pkl','rb'))
Hd=pickle.load(open('handles.pkl','rb')); valuevars=set(Hd['value'])
atomvalvars={}
for a in E.res:
    s=set()
    for u in vars_of(E.atoms[a]): s|=fa(u)
    atomvalvars[a]=[u for u in s if u in valuevars]
SHIFT=set([24468,18956])
for n in M['NODE']:
    for d in M['OUT'][n]: SHIFT.update([d['va'],d['vb'],d['vab']])
SHIFT-=set(M['live']); SHIFT-=set(M['dead'])
def relift(vv):
    r=E.run(vv); bad=[]
    for a in E.res:
        s=SL.get(a)
        if not s: continue
        i=E.residx[a]; cur=r[i]
        if cur==0: continue
        if cur%s: bad.append(a); continue
        vv[atomh[a][0]] -= cur//s
    return bad
def make(S,tag,rounds=60):
    v,isl,valn=assignment(set(S),ORIENT)
    v[24468]=T1; v[18956]=T2
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    for rd in range(rounds):
        bad=relift(vv)
        if not bad: break
        r=E.run(vv); fixed=0
        for a in bad:
            i=E.residx[a]; cur=r[i]; sm=abs(SL[a])
            if cur%p: continue                      # genuine mod-p failure, not a divisibility issue
            done=False
            imm=[x for x in vars_of(E.atoms[a]) if x in SHIFT]
            cands=imm+[x for x in atomvalvars[a] if x in SHIFT and x not in imm]
            for w in cands:
                old=vv[w]; vv[w]=old+p; d=E.run(vv)[i]-cur; vv[w]=old
                if d==0: continue
                g=gcd(d,sm)
                if cur%g: continue
                mm=sm//g
                t=(-(cur//g))*pow((d//g)%mm,-1,mm)%mm if mm>1 else 0
                vv[w]=old+p*t; done=True; fixed+=1; break
            if not done: pass
        if fixed==0: break
    bad=relift(vv)
    r=E.run(vv)
    nz=[i for i,x in enumerate(r) if x]
    print('%s |S|=%d  residual nonzero atoms %d  (unliftable %d)'%(tag,len(S),len(nz),len(bad)))
    for i in nz: print('     ',E.res[i][:140])
    out={'x_%d'%i:vv[i] for i in range(NV) if vv[i]}
    json.dump(out,open('assign_%s.json'%tag,'w'))
    print('   wrote assign_%s.json  (%d vars, max digits %d)'%(tag,len(out),max(len(str(abs(x))) for x in out.values())))
    return vv
#MAINSTART
if __name__=='__main__':
    make([M['live'][0]],'L1')
