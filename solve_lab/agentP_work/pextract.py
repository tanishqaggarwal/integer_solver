#!/usr/bin/env python3
"""Agent P: extract all 383 law-blocks, verify the template, record wiring."""
import pickle,sys,json
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
QV=24453
qpos=pickle.load(open(W+'qpos.pkl','rb'))['qpos']

def at(i): return AP[topo[i]]

def lin(ap):
    """return dict var->coef and const, or None if not linear"""
    d={};c0=0
    for m,c in ap.items():
        if len(m)==0: c0=c
        elif len(m)==1: d[m[0]]=d.get(m[0],0)+c
        else: return None,None
    return d,c0

def prod(ap):
    """ap == out - k*u*v  ->  (out,k,u,v)"""
    o=None;k=None;uv=None
    for m,c in ap.items():
        if len(m)==1 and o is None and abs(c)==1: o=(m[0],c)
        elif len(m)==2: uv=m; k=c
        else: return None
    if o is None or uv is None: return None
    # normalise: out = -k/c_o * u*v
    return (o[0], -k//o[1] if k% o[1]==0 else None, uv[0], uv[1])

def getlin(i):
    d,c0=lin(at(i))
    return d,c0

fails=[]; blocks=[]
for bi,q in enumerate(qpos):
    B={}
    try:
        # q: E = D + Q  (3 linear terms: E(+1), D(-1)... signs vary)
        d,c0=getlin(q); assert c0==0 and len(d)==3 and QV in d
        cq=d[QV]
        others=[x for x in d if x!=QV]
        # E is the one with coefficient == -cq ... both others have coefs; E defined here
        E=outof[topo[q]]; assert E in d
        Dv=[x for x in others if x!=E][0]
        assert d[E]==-cq and d[Dv]==cq, (d,E,Dv)
        B['E']=E; B['Q_off']=1 if d[E]*cq<0 else -1
        # q-1: D = C + i1
        d,_=getlin(q-1); assert len(d)==3 and Dv in d
        s=-d[Dv]
        rest=[x for x in d if x!=Dv]
        # q-2: C = i5 + i2
        Cv=None
        d2,_=getlin(q-2)
        for x in rest:
            if outof[topo[q-2]]==x: Cv=x
        assert Cv is not None
        i1=[x for x in rest if x!=Cv][0]
        assert d[i1]==s and d[Cv]==s
        assert len(d2)==3 and Cv in d2
        s2=-d2[Cv]; r2=[x for x in d2 if x!=Cv]
        # q-5: A = i1 - i2  -> identifies i2
        d5,_=getlin(q-5); A=outof[topo[q-5]]
        assert len(d5)==3 and A in d5
        sA=d5[A]; rA=[x for x in d5 if x!=A]
        assert i1 in rA, (i1,rA)
        i2=[x for x in rA if x!=i1][0]
        assert d5[i1]==-sA and d5[i2]==sA   # A = i1 - i2
        assert i2 in r2
        i5=[x for x in r2 if x!=i2][0]
        assert d2[i2]==s2 and d2[i5]==s2    # C = i2 + i5
        # q-4: B = i4 - i3
        d4,_=getlin(q-4); Bv=outof[topo[q-4]]
        assert len(d4)==3 and Bv in d4
        sB=d4[Bv]; rB=[x for x in d4 if x!=Bv]
        # q-3: A2 = A*A
        p=prod(at(q-3)); assert p and p[2]==A and p[3]==A and p[1]==1, p
        A2=p[0]
        # q+1: F = E*A2
        p=prod(at(q+1)); assert p and set((p[2],p[3]))=={E,A2} and p[1]==1
        F=p[0]
        # q+2: G = B*B
        p=prod(at(q+2)); assert p and p[2]==Bv and p[3]==Bv and p[1]==1
        G=p[0]
        # q+3: N1 = F - G
        d3,_=getlin(q+3); N1=outof[topo[q+3]]
        assert len(d3)==3 and set(d3)=={N1,F,G} and d3[F]==-d3[G]
        sN1 = -d3[F]//d3[N1]   # N1 = sN1*(F-G)
        # q+4: H = i3 + i6
        dH,_=getlin(q+4); H=outof[topo[q+4]]
        assert len(dH)==3
        rH=[x for x in dH if x!=H]
        # i3 is the one in rB
        cand=[x for x in rH if x in rB]
        assert len(cand)==1, (rH,rB)
        i3=cand[0]; i6=[x for x in rH if x!=i3][0]
        i4=[x for x in rB if x!=i3][0]
        assert d4[i4]==-sB and d4[i3]==sB   # B = i4 - i3
        assert dH[i3]==dH[i6]==-dH[H]       # H = i3 + i6
        pass
        # q+5: I = A*H
        p=prod(at(q+5)); assert p and set((p[2],p[3]))=={A,H} and p[1]==1
        I=p[0]
        # q+6: J = i2 - i5
        dJ,_=getlin(q+6); J=outof[topo[q+6]]
        assert len(dJ)==3 and dJ[i2]==-dJ[J] and dJ[i5]==dJ[J]
        # q+7: K = B*J
        p=prod(at(q+7)); assert p and set((p[2],p[3]))=={Bv,J} and p[1]==1
        K=p[0]
        # q+8: N2 = K - I
        d8,_=getlin(q+8); N2=outof[topo[q+8]]
        assert len(d8)==3 and set(d8)=={N2,K,I} and d8[K]==-d8[I]
        sN2 = -d8[K]//d8[N2]   # N2 = sN2*(K-I)
        # q+9..q+17: three outputs
        outs=[]
        for j in range(3):
            da,_=getlin(q+9+3*j); ta=outof[topo[q+9+3*j]]
            db,_=getlin(q+10+3*j); tb=outof[topo[q+10+3*j]]
            dc,_=getlin(q+11+3*j); tc=outof[topo[q+11+3*j]]
            assert set(da)=={ta,N1} and set(db)=={tb,N2}
            ca=-da[N1]//da[ta]; cb=-db[N2]//db[tb]
            assert set(dc)=={tc,ta,tb}
            sgn=dc[tc]
            outs.append((ca*(-dc[ta]//sgn), cb*(-dc[tb]//sgn), tc))
        B.update(dict(i1=i1,i2=i2,i3=i3,i4=i4,i5=i5,i6=i6,N1=N1,N2=N2,outs=outs,q=q,sN1=sN1,sN2=sN2))
        # liveness + mux: q+18..q+37
        p=prod(at(q+18)); a,b=p[2],p[3]
        B['live']=(a,b)
        mux=[]
        for base in (q+28,q+33):
            terms=[]
            for j in range(3):
                pp=prod(at(base+j)); terms.append(pp)
            mux.append(terms)
        B['mux']=mux
        B['muxout']=(outof[topo[q+32]],outof[topo[q+37]])
        blocks.append(B)
    except (AssertionError,TypeError,IndexError,KeyError) as e:
        fails.append((bi,q,repr(e)[:120]))
print("blocks matched:",len(blocks),"failed:",len(fails))
for f in fails[:10]: print("  FAIL",f)
pickle.dump(blocks,open(W+'blocks.pkl','wb'))
# coefficient matrices
mats=Counter(tuple(sorted((o[0],o[1]) for o in B['outs'])) for B in blocks)
print("distinct 3x2 output matrices:",len(mats))
qoff=Counter(B['Q_off'] for B in blocks); print('Q sign:',qoff)
print('sN1:',Counter(B['sN1'] for B in blocks));print('sN2:',Counter(B['sN2'] for B in blocks))
for s_,c in mats.most_common(5): print(c, s_)
