"""Weak-discrete-log search for k with k*G = P*  (G = bit x_2779, P* = target point)."""
import ecfast as E, ec, json, sys, time, numpy as np
p=E.p; N=E.N
P,side=ec.load(); T=ec.target()
expo=json.load(open('expo.json'))
byexp={v:int(k) for k,v in expo.items()}
G=P[byexp[0]]
PT=[P[byexp[i]] for i in range(256)]
LOG=open('dl_log.txt','a')
def say(s):
    print(s,flush=True); LOG.write(s+'\n'); LOG.flush()

def ladder(start_jac, add_affine, count, chunk=1<<16):
    """yield affine points start, start+A, start+2A, ... (count of them)"""
    J=start_jac
    buf=[]
    for i in range(count):
        buf.append(J)
        J=E.jadd_affine(J,add_affine)
        if len(buf)==chunk:
            for q in E.batch_norm(buf): yield q
            buf=[]
    if buf:
        for q in E.batch_norm(buf): yield q

def phase_runs():
    say('--- PHASE A: k = 2^a * (2^m - 1)  (runs of consecutive set bits) ---')
    t0=time.time(); hits=[]
    for a in range(256):
        S=None
        for m in range(a,256):
            S=E.aadd(S,PT[m])
            if S==T: hits.append((a,m))
    say('  runs tested %d  hits %s  %.1fs'%(256*257//2,hits,time.time()-t0))
    return hits

def phase_bsgs(mbits=22, gbits=22):
    m=1<<mbits
    say('--- PHASE B: BSGS, baby=2^%d giant=2^%d  => covers k < 2^%d ---'%(mbits,gbits,mbits+gbits))
    t0=time.time()
    keys=np.empty(m,dtype=np.int64)
    # baby[j] = j*G for j=1..m  (j=0 handled separately)
    J=(G[0],G[1],1)
    i=0
    for q in ladder(J,G,m):
        keys[i]=E.key(q[0]); i+=1
    say('  baby steps built %.1fs'%(time.time()-t0))
    order=np.argsort(keys,kind='stable'); skeys=keys[order]
    say('  sorted %.1fs'%(time.time()-t0))
    mG=E.amul(m,G); nmG=E.aneg(mG)
    found=None
    t1=time.time()
    cur=T
    # giant: T - i*m*G for i=0..2^gbits
    J=(cur[0],cur[1],1)
    gi=0
    for q in ladder(J,nmG,1<<gbits):
        kk=E.key(q[0])
        idx=np.searchsorted(skeys,kk)
        while idx<len(skeys) and skeys[idx]==kk:
            j=int(order[idx])+1
            cand=(gi*m+j)%N
            if E.amul(cand,G)==T: found=cand; break
            idx+=1
        if found is not None: break
        if q==T and gi>0: pass
        gi+=1
        if gi%(1<<19)==0: say('    giant %d/%d  %.0fs'%(gi,1<<gbits,time.time()-t1))
    say('  BSGS done %.1fs  found=%s'%(time.time()-t0,found))
    return found,keys,order

def phase_shift(keys,order):
    say('--- PHASE C: k = 2^a * s with s < 2^%d (strip trailing zeros) ---'%int(np.log2(len(keys))))
    skeys=keys[order]
    inv2=pow(2,N-2,N)
    Q=T; hits=[]
    for a in range(256):
        kk=E.key(Q[0])
        idx=np.searchsorted(skeys,kk)
        while idx<len(skeys) and skeys[idx]==kk:
            j=int(order[idx])+1
            cand=(j*pow(2,a,N))%N
            if E.amul(cand,G)==T: hits.append((a,j,cand))
            idx+=1
        Q=E.amul(inv2,Q)
    say('  shifted-target hits: %s'%hits)
    return hits

if __name__=='__main__':
    what=sys.argv[1] if len(sys.argv)>1 else 'all'
    if what in ('all','A'): phase_runs()
    if what in ('all','B','BC'):
        mb=int(sys.argv[2]) if len(sys.argv)>2 else 22
        gb=int(sys.argv[3]) if len(sys.argv)>3 else 22
        f,keys,order=phase_bsgs(mb,gb)
        if f is not None: say('*** DISCRETE LOG FOUND: k = %d ***'%f)
        else: phase_shift(keys,order)
