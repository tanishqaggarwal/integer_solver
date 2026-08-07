#!/usr/bin/env python3
"""Agent P: independent fold evaluator over the recovered 383-stage network."""
import pickle,sys,json
from collections import Counter,defaultdict
sys.set_int_max_str_digits(10**7)
W='/home/user/integer_solver/solve_lab/agentP_work/'
P=115792089237316195423570985008687907853269984665640564039457584007908834671663
Q=97553848499418123410591666447050222001188385549510401465815187079080512838891
D=pickle.load(open(W+'model4.pkl','rb')); AP=D['AP']
S=pickle.load(open(W+'slp.pkl','rb')); topo=S['topo']; outof=S['outof']
G=pickle.load(open(W+'graph.pkl','rb')); stages=G['stages']; outmap=G['outmap']
Wd=pickle.load(open(W+'wire3.pkl','rb')); par=Wd['par']; leafcoord=Wd['leafcoord']; Z=set(Wd['Z'])
T=pickle.load(open(W+'topo.pkl','rb')); src=T['src']; supp=T['supp']
LEAVES=pickle.load(open(W+'leaves.pkl','rb'))
R=lambda x:par[x]
SEL=sorted(set(a for a,b,k in LEAVES))
SIDX={s:i for i,s in enumerate(SEL)}
LK=defaultdict(list)
for a,b,k in LEAVES: LK[a].append((b,k%P))
# order the two coords of each leaf by the (i2,i3) / (i1,i4) role at its consuming stage
def leafpair(selvar, stage_j, slot):
    s=stages[stage_j]; pr=s['X'] if slot==0 else s['Y']
    out=[]
    for c in pr:
        for (b,k) in LK[selvar]:
            if R(b)==R(c): out.append(k)
    return tuple(out)

# resolve sources: for each stage, (kindX, refX), (kindY, refY)
SRC=[]
for j,row in enumerate(src):
    r=[]
    for slot,(k,v) in enumerate(row):
        if k=='S': r.append(('S',v))
        elif k=='L': r.append(('L',v,leafpair(v,j,slot)))
        else: r.append(('0',None))
    SRC.append(r)
print("resolved:",Counter(x[0] for r in SRC for x in r))

# root: block index 382 (q=19004) uses stage 381 & 380 as X/Y
ROOT_X=381; ROOT_Y=380
TGT_X=91416258160755509149180373473728639746431157665678710450404458852172057265575180278101002 % P
TGT_Y=1114942656963403660822546820446916783439088877768247923308647546252105232931473698035897478439338 * pow(8863713,P-2,P) % P

def law(X,Y):
    (x1,y1),(x2,y2)=X,Y
    u=(x2-x1)%P; v=(y2-y1)%P
    if u==0: return None            # degenerate
    lam=v*pow(u,P-2,P)%P
    x3=(lam*lam-x1-x2-Q)%P
    y3=(lam*(x1-x3)-y1)%P
    return (x3,y3)

def fold(bits, root_only=True):
    """bits: list of 256 ints. returns (pair,live) at root."""
    memo={}
    def ev(j):
        if j in memo: return memo[j]
        vals=[]
        for k in SRC[j]:
            if k[0]=='S': vals.append(ev(k[1]))
            elif k[0]=='L':
                b=bits[SIDX[k[1]]]
                vals.append((k[2],1) if b else ((0,0),0))
            else: vals.append(((0,0),0))
        (X,a),(Y,b)=vals
        if a and b:
            r=law(X,Y)
            res=(r,1) if r else (None,1)
        elif a: res=(X,1)
        elif b: res=(Y,1)
        else: res=((0,0),0)
        memo[j]=res; return res
    sys.setrecursionlimit(10000)
    X,a=ev(ROOT_X); Y,b=ev(ROOT_Y)
    if a and b:
        if X is None or Y is None: return None,1
        return law(X,Y),1
    if a: return X,1
    if b: return Y,1
    return (0,0),0

if __name__=='__main__':
    g=[0]*38748
    for k,v in json.load(open('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')).items(): g[int(k[2:])]=int(v)
    bits=[1 if g[s] else 0 for s in SEL]
    on=[i for i,b in enumerate(bits) if b]
    print("deliverable ON leaves:",on)
    r,live=fold(bits)
    print("fold root =",r)
    print("deliverable root vars: x13682 =",g[24468]%P if g[24468] else g[13682]%P)
    # find actual values of root class vars
    print("  g[x13682]%P =",g[13682]%P, " g[x37892]%P =",g[37892]%P)
    print("TARGET =",(TGT_X,TGT_Y))
    print("MATCH root==deliverable?", r==(g[13682]%P,g[37892]%P))
    print("MATCH root==target?", r==(TGT_X,TGT_Y))
