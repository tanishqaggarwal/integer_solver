"""MITM for low Hamming weight: all subsets of size<=3 of the 256 exponents, matched against
   T - (subsets of size<=3).  Covers weight <= 6 exactly."""
import ecfast as E, ec, json, sys, time, numpy as np
p=E.p; N=E.N
P,side=ec.load(); T=ec.target()
expo=json.load(open('expo.json')); byexp={v:int(k) for k,v in expo.items()}
PT=[P[byexp[i]] for i in range(256)]
LOG=open('mitm_log.txt','a')
def say(s):
    print(s,flush=True); LOG.write(s+'\n'); LOG.flush()

def build_sums(maxw=3):
    """returns (keys np.int64, meta list of tuples) for every subset of size 1..maxw."""
    t0=time.time()
    jac=[]; meta=[]
    # size 1
    for i in range(256):
        jac.append((PT[i][0],PT[i][1],1)); meta.append((i,))
    # size 2 (affine, only 32640)
    S2=[]
    for i in range(256):
        for j in range(i+1,256):
            S2.append(((i,j)))
    # compute pair sums in jacobian from P_i + P_j
    pj=[]
    for (i,j) in S2:
        pj.append(E.jadd_affine((PT[i][0],PT[i][1],1),PT[j]))
    pa=E.batch_norm(pj)
    for idx,(i,j) in enumerate(S2):
        jac.append((pa[idx][0],pa[idx][1],1)); meta.append((i,j))
    say('  size<=2 built %d  %.1fs'%(len(meta),time.time()-t0))
    if maxw>=3:
        buf=[]; bmeta=[]
        keys3=[]; meta3=[]
        CH=1<<16
        for idx,(i,j) in enumerate(S2):
            base=pa[idx]
            bj=(base[0],base[1],1)
            for k in range(j+1,256):
                buf.append(E.jadd_affine(bj,PT[k])); bmeta.append((i,j,k))
                if len(buf)==CH:
                    for q,mm in zip(E.batch_norm(buf),bmeta):
                        jac.append((q[0],q[1],1)); meta.append(mm)
                    buf=[]; bmeta=[]
        if buf:
            for q,mm in zip(E.batch_norm(buf),bmeta):
                jac.append((q[0],q[1],1)); meta.append(mm)
    say('  size<=%d built %d subsets  %.1fs'%(maxw,len(meta),time.time()-t0))
    keys=np.fromiter((E.key(t[0]) for t in jac),dtype=np.int64,count=len(jac))
    xs=[t[0] for t in jac]; ys=[t[1] for t in jac]
    return keys,meta,xs,ys

if __name__=='__main__':
    t0=time.time()
    keys,meta,xs,ys=build_sums(3)
    order=np.argsort(keys,kind='stable'); skeys=keys[order]
    say('sorted %d  %.1fs'%(len(keys),time.time()-t0))
    # weight<=3 direct
    hit=[]
    kt=E.key(T[0])
    idx=np.searchsorted(skeys,kt)
    while idx<len(skeys) and skeys[idx]==kt:
        i=int(order[idx])
        if (xs[i],ys[i])==T: hit.append(meta[i])
        idx+=1
    say('weight<=3 hits: %s'%hit)
    # weight<=6: for each subset S, look up T - sum(S)
    found=[]
    CH=1<<16
    for start in range(0,len(meta),CH):
        end=min(start+CH,len(meta))
        js=[E.jadd_affine((T[0],T[1],1),(xs[i],(-ys[i])%p)) for i in range(start,end)]
        aff=E.batch_norm(js)
        for off,q in enumerate(aff):
            if q is None: continue
            kk=E.key(q[0])
            ii=np.searchsorted(skeys,kk)
            while ii<len(skeys) and skeys[ii]==kk:
                j=int(order[ii])
                if (xs[j],ys[j])==q:
                    a=set(meta[start+off]); b=set(meta[j])
                    if not (a&b): found.append((sorted(a),sorted(b)))
                ii+=1
        if start % (1<<20) == 0: say('  scanned %d/%d  %.0fs  found=%d'%(start,len(meta),time.time()-t0,len(found)))
        if found: break
    say('weight<=6 hits: %s   total %.1fs'%(found[:3],time.time()-t0))
