import sys, pickle, time, collections
src=open('/home/user/integer_solver/solve_lab/agentT_work/mirror/L/lift.py').read().split("if __name__")[0]
exec(src)
t0=time.time()
vv=[0]*NV; base=E.run(vv)
SL={}
for k,a in enumerate(E.res):
    hs=atomh[a]
    if len(hs)!=1: continue
    h=hs[0]; i=E.residx[a]
    vv[h]=1; s=E.run(vv)[i]-base[i]; vv[h]=0
    SL[a]=s
    if k%1000==0: print(k,time.time()-t0,flush=True)
print('slopes for %d atoms, %.0fs'%(len(SL),time.time()-t0))
print('slope/p hist', collections.Counter(abs(s)//p for s in SL.values() if s%p==0).most_common(6))
print('slopes not divisible by p:',sum(1 for s in SL.values() if s%p))
pickle.dump(SL,open('slopes.pkl','wb'))
