import sys, json, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from close4 import *
P=2**256-2**32-977
K1=125787314747601108116039725163361763116550465675981151838811516827327919228823597744635626
K2=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002
W='/home/user/integer_solver/solve_lab/agentC_work/'
BI=json.load(open(W+'bitinfo.json'))
def plan(target_bit, partner):
    d=BI[str(target_bit)]; q=BI[str(partner)]
    Tx=int(q['px']); Ty=int(q['py'])
    ctrl={22162:K2,30213:K1,target_bit:1,partner:1}
    det=[]
    for pin,T in ((d['xpin'],Tx),(d['ypin'],Ty)):
        X=pin['X']; C=pin['C']; m=pin['mult']; H=pin['H']
        Dv=T-C
        # choose t with m | (D + p*t)
        if m==1: t=0
        else:
            g=pow(P%m,-1,m) if __import__('math').gcd(P%m,m)==1 else None
            if g is None: return None
            t=(-Dv*g)%m
        Xv=T+P*t
        assert (Xv-C)%m==0
        ctrl[X]=Xv; ctrl[H]=(Xv-C)//m
        det.append(H)
    return ctrl,det
if __name__=='__main__':
    tb=int(sys.argv[1]); partners=[int(x) for x in sys.argv[2:]]
    best=0
    for pa in partners:
        r=plan(tb,pa)
        if r is None: print('skip',pa); continue
        ctrl,det=r
        t=time.time()
        sc,v,nz=closure4(ctrl,detach=det,rounds=12)
        print('target x_%d partner x_%d -> score %d  nz=%s  %.1fs'%(tb,pa,sc,[(a,len(L.atom2eq.get(a,{}))) for a in nz],time.time()-t),flush=True)
        if sc>best:
            best=sc
            json.dump({f'x_{i}':v[i] for i in range(L.NVARS) if v[i]!=0},open(W+'P10513_%d_%d.json'%(sc,pa),'w'))
    print('BEST',best)
