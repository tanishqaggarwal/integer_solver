import sys, pickle, json, time, collections
src=open('/home/user/integer_solver/solve_lab/agentL_work/lift.py').read().split("if __name__")[0]
exec(src)
SL=pickle.load(open('slopes.pkl','rb'))
def make(S,tag):
    v,isl,valn=assignment(set(S),ORIENT)
    v[24468]=T1; v[18956]=T2
    vv=[0]*NV
    for k,x in v.items(): vv[k]=x
    r=E.run(vv)
    bad=[]
    for a in E.res:
        s=SL.get(a)
        if s is None or s==0: continue
        i=E.residx[a]; cur=r[i]
        if cur==0: continue
        if cur%s: bad.append(a); continue
        vv[atomh[a][0]] -= cur//s
    r=E.run(vv)
    nz=[i for i,x in enumerate(r) if x]
    nzp=[i for i,x in enumerate(r) if x%p]
    print('%s |S|=%d  unliftable %d  nonzero atoms %d (mod p: %d)  root=%s'%(tag,len(S),len(bad),len(nz),len(nzp),valn[ROOT]))
    for i in nz: print('     ',E.res[i][:140])
    out={'x_%d'%i:vv[i] for i in range(NV) if vv[i]}
    json.dump(out,open('assign_%s.json'%tag,'w'))
    return vv,r
if __name__=='__main__':
    make([M['live'][0]],'L1')
    make(M['live'][:2],'L2')
