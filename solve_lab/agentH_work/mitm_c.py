"""Complement MITM: find S with sum_{i in S} P_i = ALL - T  (i.e. k of Hamming weight >= 250)."""
import ecfast as E, ec, json, time, numpy as np
import mitm as M   # reuses build_sums; but override target
p=E.p
PT=M.PT
ALL=None
for i in range(256): ALL=E.aadd(ALL,PT[i])
T2=E.aadd(ALL,E.aneg(M.T))
print('complement target computed',flush=True)
t0=time.time()
keys,meta,xs,ys=M.build_sums(3)
order=np.argsort(keys,kind='stable'); skeys=keys[order]
print('sorted %d %.1fs'%(len(keys),time.time()-t0),flush=True)
hit=[]
kt=E.key(T2[0]); idx=np.searchsorted(skeys,kt)
while idx<len(skeys) and skeys[idx]==kt:
    i=int(order[idx])
    if (xs[i],ys[i])==T2: hit.append(meta[i])
    idx+=1
print('complement weight<=3 hits:',hit,flush=True)
found=[]; CH=1<<16
for start in range(0,len(meta),CH):
    end=min(start+CH,len(meta))
    js=[E.jadd_affine((T2[0],T2[1],1),(xs[i],(-ys[i])%p)) for i in range(start,end)]
    for off,q in enumerate(E.batch_norm(js)):
        if q is None: continue
        kk=E.key(q[0]); ii=np.searchsorted(skeys,kk)
        while ii<len(skeys) and skeys[ii]==kk:
            j=int(order[ii])
            if (xs[j],ys[j])==q:
                a=set(meta[start+off]); b=set(meta[j])
                if not (a&b): found.append((sorted(a),sorted(b)))
            ii+=1
    if found: break
print('complement weight<=6 hits:',found[:3],'  %.1fs'%(time.time()-t0),flush=True)
