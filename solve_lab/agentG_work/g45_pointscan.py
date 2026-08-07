"""Per-bit: solved six coordinates -> secp256k1 points; tests F_p-affinity of (A,B)
and GROUP-affinity of D = P3 - (P1+P2)."""
import os, sys, itertools, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F, gpt
import gsym2 as G
from gsym2 import L, ad, P
NB=F.NB; n=len(NB); ixm={u:i for i,u in enumerate(NB)}

def frame(FL):
    v=list(F.v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    r=F.analyse(v)
    val,_=G.build(v,NB,cap=6)
    f1=G.evalatom(19297,val,6); f2=G.evalatom(19299,val,6)
    o={'flips':FL,'ninc':len(r['incchecks']),'nzc':len(r['nzc']),'nres':len(r['res'])}
    if isinstance(f1,int) or isinstance(f2,int):
        o['core']='DEAD'; return o
    Ap,Bp=gpt.pencil(f1,f2); lab=gpt.label(Ap,Bp,NB); o['label']=lab
    if lab is None: o['core']='UNLABELLED'; return o
    piv,R=r['piv'],r['R']
    t={c:(v[NB[c]]%P) for c in range(n) if c not in piv}
    def value(u):
        j=ixm.get(u)
        if j is None: return v[u]%P
        if j not in piv: return t[j]
        row=R[piv[j]]; val_=row.get(n,0)%P
        for c,vv in row.items():
            if c!=n and c!=j: val_=(val_-vv*t.get(c,0))%P
        return val_
    co={k:value(u) for k,u in lab.items()}
    o['coord']=co
    o['free']={k:(ixm.get(u) not in piv) for k,u in lab.items()}
    P1=gpt.tosec(co['x1'],co['y1']); P2=gpt.tosec(co['x2'],co['y2']); P3=gpt.tosec(co['x3'],co['y3'])
    o['pts']=(P1,P2,P3)
    o['oncurve']=[(q[1]*q[1]-pow(q[0],3,P)-7)%P==0 for q in (P1,P2,P3)]
    o['D']=gpt.sub(P3,gpt.add(P1,P2)) if all(o['oncurve']) else None
    # A,B evaluated at the solved point
    def ev(f):
        s=0
        for m,c in f.items():
            tt=c
            for k,e in m: tt=tt*pow(value(NB[k]),e,P)%P
            s=(s+tt)%P
        return s
    o['AB']=(ev(Ap),ev(Bp))
    return o

if __name__=='__main__':
    BITS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [47,91,112,438,490,542,853,1203,1357,1413]
    base=frame([]); res={():base}
    print('BASE oncurve=%s AB=%s'%(base.get('oncurve'),str(base['AB'])[:50]))
    print('  D =',base.get('D'),flush=True)
    for b in BITS:
        o=frame([b]); res[(b,)]=o
        print('x%-6d nres=%-2d label=%s oncurve=%s free=%s\n        D=%s'%(
            b,o['nres'],{k:o.get('label',{}).get(k) for k in ['x1','y1','x2','y2','x3','y3']} if o.get('label') else None,
            o.get('oncurve'),[k for k,f in o.get('free',{}).items() if f],str(o.get('D'))[:80]),flush=True)
    print('\n--- pair tests (F_p-affine? group-affine?) ---')
    for i,j in itertools.combinations(BITS[:6],2):
        o=frame([i,j]); res[(i,j)]=o
        li,lj=res[(i,)],res[(j,)]; msg=[]
        if o.get('AB') and li.get('AB') and lj.get('AB'):
            pred=((li['AB'][0]+lj['AB'][0]-base['AB'][0])%P,(li['AB'][1]+lj['AB'][1]-base['AB'][1])%P)
            msg.append('Fp-AFFINE' if pred==o['AB'] else 'Fp-nonlinear')
        if o.get('D') and li.get('D') and lj.get('D') and base.get('D'):
            predD=gpt.sub(gpt.add(li['D'],lj['D']),base['D'])
            msg.append('GROUP-AFFINE' if predD==o['D'] else 'group-nonlinear')
        print('  (%d,%d): nres=%d oncurve=%s  %s'%(i,j,o['nres'],o.get('oncurve'),' '.join(msg)),flush=True)
    pickle.dump(res,open('/home/user/integer_solver/solve_lab/agentG_work/pointscan.pkl','wb'))
