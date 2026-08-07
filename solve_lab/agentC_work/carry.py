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
E1=[i for i,b in enumerate(chain) if side[chain[i]]=='s1']
E2=[i for i,b in enumerate(chain) if side[chain[i]]=='s2']
print('|E1|',len(E1),'|E2|',len(E2),'union',len(set(E1)|set(E2)))
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
S1=set(E1); S2=set(E2)
def solve(n):
    m=0; c=0; A=[]; B=[]
    for e in range(256):
        nb=(n>>e)&1
        if e in S2:
            me=(nb+c)&1
            if me: B.append(e); m|=1<<e
            s=nb+me+c
            r=s&1; c=s>>1
            assert r==0
        else:
            s=nb+c; r=s&1; c=s>>1
            if r: A.append(e)
    return c,A,B,m
for sgn,n in [(+1,N)]:
    c,A,B,m=solve(n)
    kA=sum(1<<e for e in A); kB=sum(1<<e for e in B)
    print('final carry',c,'|A|',len(A),'|B|',len(B))
    print('kA-kB == N ?', kA-kB==n, ' (kA<2^256:',kA<2**256,')')
    if c==0 and kA-kB==n:
        json.dump({'A':A,'B':B},open('/home/user/integer_solver/solve_lab/agentC_work/AB.json','w'))
        # verify on the curve
        G=pts[chain[0]]
        PA=None
        for e in A: PA=add(PA,pts[chain[e]])
        PB=None
        for e in B: PB=add(PB,pts[chain[e]])
        print('P_A == P_B ?', PA==PB)
        print('PA',PA)
