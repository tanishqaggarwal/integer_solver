"""U5: independent recovery of the curve, the 256 leaf points, and the leaf->exponent map.
Uses ONLY my own parse (v_defs.pkl).  No other agent's data."""
import pickle, collections
B='/home/user/integer_solver/solve_lab/agentU_work/'
D=pickle.load(open(B+'v_defs.pkl','rb'))
CONST=D['CONST']; LEAFPIN=D['LEAFPIN']
p=CONST[26064]
print('p bitlen',p.bit_length(),'prime?',pow(2,p-1,p)==1)
# group leaf pins by selector
bysel=collections.defaultdict(list)
for sel,w,C,z,m in LEAFPIN: bysel[sel].append((w,C%p))
assert all(len(v)==2 for v in bysel.values()), 'not 2 pins/sel'
sels=sorted(bysel)
print('selectors',len(sels))
# ---- solve for the shift s and b algebraically: y^2 = (x+s)^3 + b ----
# orientation (which of a selector's two pin constants is X) is unknown per selector.
import itertools
def quad(P1,P2):
    (x1,y1),(x2,y2)=P1,P2
    A=3*(x1-x2)%p; Bq=3*(x1*x1-x2*x2)%p
    C=((x1*x1*x1-x2*x2*x2)-(y1*y1-y2*y2))%p
    return A,Bq,C
c3=[[bysel[s][0][1],bysel[s][1][1]] for s in sels[:3]]
found=None
for o in itertools.product((0,1),repeat=3):
    P=[(c3[i][o[i]], c3[i][1-o[i]]) for i in range(3)]
    A1,B1,C1=quad(P[0],P[1]); A2,B2,C2=quad(P[0],P[2])
    den=(A2*B1-A1*B2)%p
    if den==0: continue
    s=(-(A2*C1-A1*C2))*pow(den,p-2,p)%p
    b=(P[0][1]**2-pow((P[0][0]+s)%p,3,p))%p
    ok=True; assign={}
    for sl in sels:
        c1,c2=bysel[sl][0][1],bysel[sl][1][1]
        hit=[k for k,(X,Y) in enumerate([(c1,c2),(c2,c1)]) if (Y*Y-pow((X+s)%p,3,p)-b)%p==0]
        if not hit: ok=False; break
        assign[sl]=hit[0] if len(hit)==1 else hit[0]
        if len(hit)>1: print('ambiguous orientation at sel',sl)
    if ok: found=(s,b,assign); break
assert found, 'NOT FOUND'
shift,b,ORI=found
print('shift =',shift); print('b =',b); print('3*shift mod p =',3*shift%p)
pts={s:((bysel[s][ORI[s]][1]+shift)%p, bysel[s][1-ORI[s]][1]) for s in sels}
print('all 256 on Y^2=X^3+b :', all((Y*Y-pow(X,3,p)-b)%p==0 for X,Y in pts.values()))
# ---- doubling chain, in the shifted (a=0) Weierstrass frame ----
def dbl(P):
    X,Y=P
    if Y%p==0: return None
    l=3*X*X%p*pow(2*Y,p-2,p)%p
    x3=(l*l-2*X)%p; return (x3,(l*(X-x3)-Y)%p)
inv={v:k for k,v in pts.items()}
assert len(inv)==256
nxt={}; hits=0
for s,P in pts.items():
    d=dbl(P)
    if d in inv: nxt[s]=inv[d]; hits+=1
print('doublings landing inside the leaf set: %d/256'%hits)
src=set(pts)-set(nxt.values())
print('leaves that are nobody\'s double:',len(src), 'leaves with no double in set:', len(set(pts)-set(nxt)))
base=list(src)[0]
exp={}; cur=base; e=0
while True:
    exp[cur]=e
    if cur not in nxt: break
    cur=nxt[cur]; e+=1
print('chain length',len(exp),'max exponent',max(exp.values()))
assert len(exp)==256
# group order N: verify by scalar multiplication of the base
def add(P,Q):
    if P is None: return Q
    if Q is None: return P
    X1,Y1=P; X2,Y2=Q
    if X1==X2:
        if (Y1+Y2)%p==0: return None
        return dbl(P)
    l=(Y2-Y1)*pow(X2-X1,p-2,p)%p
    x3=(l*l-X1-X2)%p; return (x3,(l*(X1-x3)-Y1)%p)
def mul(k,P):
    R=None; Q=P
    while k:
        if k&1: R=add(R,Q)
        Q=dbl(Q) if Q else None; k>>=1
    return R
N=115792089237316195423570985008687907852837564279074904382605163141518161494337
G=pts[base]
print('N*G == O :', mul(N,G) is None)
print('leaf(e) == 2^e*G for all e :', all(mul(1<<exp[s],G)==pts[s] for s in sels))
sel2exp={s:exp[s] for s in sels}
pickle.dump({'p':p,'K':3*shift%p,'shift':shift,'b':b,'ORI':ORI,'pts':pts,'sel2exp':sel2exp,
             'bysel':{s:[w for w,_ in bysel[s]] for s in sels},'N':N}, open(B+'v_leaves.pkl','wb'))
print('base selector x%d'%base)
