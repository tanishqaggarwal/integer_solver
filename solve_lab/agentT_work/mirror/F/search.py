#!/usr/bin/env python3
"""Min-conflicts / annealing search on the full atom model."""
import sys,os,json,random,time,math,collections
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
from full import Full
from fwd import NV

def roots_of(F,e,i):
    """S_e as a function of t=v[i]: fit quadratic exactly, return integer roots."""
    old=F.v[i]
    ys=[]
    for dt in (0,1,2):
        F.setvar(i,old+dt); ys.append(F.ev[e])
    F.setvar(i,old)
    c0=ys[0]; d1=ys[1]-ys[0]; d2=ys[2]-2*ys[1]+ys[0]
    a=d2//2 if d2%2==0 else None
    out=[]
    if d2==0:
        if d1==0: return []
        if c0 % d1==0: out.append(old - c0//d1)
        return out
    if a is None: return []
    b=d1-a
    # a*t^2 + b*t + c0 = 0 with t = delta
    disc=b*b-4*a*c0
    if disc<0: return []
    s=math.isqrt(disc)
    if s*s!=disc: return []
    for sg in (1,-1):
        num=-b+sg*s
        if num % (2*a)==0: out.append(old+num//(2*a))
    return out

def search(F, iters=200000, seed=0, log=None, best_cb=None):
    rnd=random.Random(seed)
    best=F.nfail; bestv=list(F.v)
    t0=time.time(); last=t0
    stall=0
    for it in range(iters):
        fails=F.fails()
        if not fails: return 0,list(F.v)
        e=rnd.choice(fails)
        vars_=set()
        for k,j in F.eqrows[e]: vars_.update(F.avars[j])
        vars_=list(vars_); rnd.shuffle(vars_)
        cand=[]
        for i in vars_[:40]:
            for t in roots_of(F,e,i):
                old=F.v[i]
                nf=F.setvar(i,t)
                cand.append((nf,i,t))
                F.setvar(i,old)
        if not cand:
            i=rnd.choice(vars_); F.setvar(i,F.v[i]+rnd.choice([-1,1])); continue
        cand.sort()
        nf,i,t=cand[0]
        if nf<=F.nfail or rnd.random()<0.05:
            F.setvar(i,t)
        else:
            nf2,i2,t2=rnd.choice(cand); F.setvar(i2,t2)
        if F.nfail<best:
            best=F.nfail; bestv=list(F.v); stall=0
            if best_cb: best_cb(best,bestv)
            if log: print('it',it,'NEW BEST nfail',best,'score',39033-best,'t',round(time.time()-t0,1),flush=True)
        else: stall+=1
        if log and time.time()-last>60:
            last=time.time(); print('it',it,'cur',F.nfail,'best',best,'t',round(time.time()-t0,1),flush=True)
        if stall>4000:
            # restart from best with random kick
            for i in range(NV): F.v[i]=bestv[i]
            F.init(F.v)
            for _ in range(3):
                fs=F.fails(); 
                if not fs: break
                e2=rnd.choice(fs); vs=set()
                for k,j in F.eqrows[e2]: vs.update(F.avars[j])
                i2=rnd.choice(list(vs)); F.setvar(i2,F.v[i2]+rnd.randrange(-3,4))
            stall=0
    return best,bestv

if __name__=='__main__':
    src=sys.argv[1] if len(sys.argv)>1 else os.path.join(HERE,'..','best','new_instance_partial_39026.json')
    seed=int(sys.argv[2]) if len(sys.argv)>2 else 0
    tag=sys.argv[3] if len(sys.argv)>3 else 'A'
    F=Full()
    d=json.load(open(src)); v=[0]*NV
    for k,x in d.items(): v[int(k[2:]) if k.startswith('x_') else int(k)]=int(x)
    print('start nfail',F.init(v),flush=True)
    BEST=[F.nfail]
    def cb(b,bv):
        if b<BEST[0]:
            BEST[0]=b
            json.dump({'x_%d'%i:bv[i] for i in range(NV) if bv[i]},open(os.path.join(HERE,'SEARCH_%s_%d.json'%(tag,39033-b)),'w'))
    b,bv=search(F,iters=10**9,seed=seed,log=True,best_cb=cb)
    print('final best',b,'score',39033-b)
