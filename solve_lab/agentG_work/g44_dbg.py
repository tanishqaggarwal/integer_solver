import os, sys
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import g29_frame as F, gpt
import gsym2 as G
from gsym2 import L, ad, P
NB=F.NB
for arg in sys.argv[1:]:
    FL=[int(x) for x in arg.split(',') if x] if arg!='-' else []
    v=list(F.v0)
    for b in FL: v[b]=1-v[b]
    ad.fwd(v,rounds=8)
    val,_=G.build(v,NB,cap=6)
    f1=G.evalatom(19297,val,6); f2=G.evalatom(19299,val,6)
    if isinstance(f1,int) or isinstance(f2,int):
        print(arg,'a19297/a19299 CONSTANT'); continue
    Ap,Bp=gpt.pencil(f1,f2)
    def ms(m): return '*'.join('x%d'%NB[k]+('^%d'%e if e>1 else '') for k,e in m) or '1'
    def sc(c):
        c%=P; return str(c) if c<10**10 else ('(-%d)'%(P-c) if P-c<10**10 else 'C')
    print('=== flips',FL)
    print(' A:',' '.join('%s*%s'%(sc(c),ms(m)) for m,c in sorted(Ap.items(),key=lambda kv:-G.mdeg(kv[0]))))
    print(' B:',' '.join('%s*%s'%(sc(c),ms(m)) for m,c in sorted(Bp.items(),key=lambda kv:-G.mdeg(kv[0]))))
    print(' label:',gpt.label(Ap,Bp,NB))
