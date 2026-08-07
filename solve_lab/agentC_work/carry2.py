import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentC_work')
from ec import *
from ort import *
chain=[int(b) for b in json.load(open('/home/user/integer_solver/solve_lab/agentC_work/chain.json'))]
pts=leafpoints()
side={}
for r,tag in [(8599,'s1'),(21839,'s1'),(25956,'s2'),(7304,'s2')]:
    for x in leaves(r):
        if x in pts: side[x]=tag
S1={i for i in range(256) if side[chain[i]]=='s1'}
S2={i for i in range(256) if side[chain[i]]=='s2'}
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
def chase(n, ctrl, target):
    """n + sum_{ctrl bits} 2^e = result, with result bits forced 0 at positions NOT in target.
       ctrl = positions where we may add a bit. Returns (carry, chosen, resultbits)."""
    ch=[]; res=[]; c=0; bad=[]
    for e in range(256):
        nb=(n>>e)&1
        if e in ctrl:
            me=(nb+c)&1          # force result bit 0
            if me: ch.append(e)
            s=nb+me+c; r=s&1; c=s>>1
            assert r==0
        else:
            s=nb+c; r=s&1; c=s>>1
            if r:
                if e in target: res.append(e)
                else: bad.append(e)
    return c,ch,res,bad
# case w=+N : N + kB = kA, ctrl=S2 (choose B), result support must be in S1
c,B,A,bad=chase(N,S2,S1)
print('w=+N : carry',c,'bad',len(bad),'|A|',len(A),'|B|',len(B))
# case w=-N : N + kA = kB, ctrl=S1 (choose A), result support must be in S2
c2,A2,B2,bad2=chase(N,S1,S2)
print('w=-N : carry',c2,'bad',len(bad2),'|A|',len(A2),'|B|',len(B2))
for tag,(cc,AA,BB,bb) in [('+N',(c,A,B,bad)),('-N',(c2,A2,B2,bad2))]:
    if cc==0 and not bb:
        kA=sum(1<<e for e in AA); kB=sum(1<<e for e in BB)
        print(tag,'VALID: kA-kB =',kA-kB, 'N=',N, kA-kB==N or kB-kA==N)
        PA=None
        for e in AA: PA=add(PA,pts[chain[e]])
        PB=None
        for e in BB: PB=add(PB,pts[chain[e]])
        print('   P_A == P_B ?', PA==PB, 'nonempty', len(AA)>0, len(BB)>0)
        json.dump({'A':sorted(AA),'B':sorted(BB),'case':tag},open('/home/user/integer_solver/solve_lab/agentC_work/AB.json','w'))
