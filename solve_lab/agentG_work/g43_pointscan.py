"""Per-bit: the exact solved six coordinates -> secp256k1 points, and D = P3-(P1+P2).
Tests (a) F_p-affinity of (A,B) and (b) GROUP-affinity of D."""
import os, sys, itertools, pickle, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F, gpt
import gsym2 as G
from gsym2 import L, ad, P
NB=F.NB; n=len(NB); ix={u:i for i,u in enumerate(NB)}

def frame(FL):
    v=list(F.v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    r=F.analyse(v)
    out={'flips':FL,'ninc':len(r['incchecks']),'nzc':len(r['nzc']),'nres':len(r['res'])}
    d=dict(r['res'])
    a1=d.get(19297,0); a2=d.get(19299,0)
    if isinstance(a1,int) and isinstance(a2,int):
        A=(gpt.c22*a1-gpt.c12*a2)*gpt._inv%P
        B=(-gpt.c21*a1+gpt.c11*a2)*gpt._inv%P
        out['AB']=(A,B)
    else: out['AB']=None
    # solved coordinate values: from the linear solution with free params at 0
    val,_=G.build(v,NB,cap=6)
    f1=G.evalatom(19297,val,6); f2=G.evalatom(19299,val,6)
    if isinstance(f1,int) or isinstance(f2,int): out['pts']=None; return out,r
    Ap,Bp=gpt.pencil(f1,f2)
    lab=gpt.label(Ap,Bp,NB)
    out['label']=lab
    if lab is None: out['pts']=None; return out,r
    piv,R=r['piv'],r['R']
    vals={}
    for k,u in lab.items():
        j=ix.get(u)
        if j is None: vals[k]=None; continue
        if j in piv:
            row=R[piv[j]]
            vals[k]=row.get(n,0)%P
            vals[k+'_free']=[c for c in row if c!=n and c!=j]
        else:
            vals[k]='FREE'
    out['coord']=vals
    try:
        P1=gpt.tosec(vals['x1'],vals['y1']); P2=gpt.tosec(vals['x2'],vals['y2']); P3=gpt.tosec(vals['x3'],vals['y3'])
        out['pts']=(P1,P2,P3)
        out['oncurve']=[(q[1]*q[1]-pow(q[0],3,P)-7)%P==0 for q in (P1,P2,P3)]
        out['D']=gpt.sub(P3,gpt.add(P1,P2))
    except Exception as e:
        out['pts']=None
    return out,r

if __name__=='__main__':
    BITS=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [47,91,112,438,490,542,853,1203]
    base,_=frame([])
    print('BASE AB=%s'%(str(base['AB'])[:60]))
    print('  label',base.get('label'))
    print('  oncurve',base.get('oncurve'),' D=',str(base.get('D'))[:70],flush=True)
    res={():base}
    for b in BITS:
        o,_=frame([b]); res[(b,)]=o
        print('bit x%-6d ninc=%d nzc=%d nres=%d oncurve=%s AB=%s D=%s'%(
            b,o['ninc'],o['nzc'],o['nres'],o.get('oncurve'),str(o['AB'])[:44],str(o.get('D'))[:44]),flush=True)
    print('\n--- pair tests ---')
    for i,j in itertools.combinations(BITS[:5],2):
        o,_=frame([i,j]); res[(i,j)]=o
        li,lj=res[(i,)],res[(j,)]
        msg=[]
        if o['AB'] and li['AB'] and lj['AB'] and base['AB']:
            pred=((li['AB'][0]+lj['AB'][0]-base['AB'][0])%P,(li['AB'][1]+lj['AB'][1]-base['AB'][1])%P)
            msg.append('Fp-affine' if pred==o['AB'] else 'Fp-NONLINEAR')
        if o.get('D') is not None and li.get('D') is not None and lj.get('D') is not None:
            predD=gpt.sub(gpt.add(li['D'],lj['D']),base['D'])
            msg.append('group-affine' if predD==o['D'] else 'group-NONLINEAR')
        print('  (%d,%d): ninc=%d nres=%d  %s'%(i,j,o['ninc'],o['nres'],' '.join(msg)),flush=True)
    pickle.dump(res,open('/home/user/integer_solver/solve_lab/agentG_work/pointscan.pkl','wb'))
