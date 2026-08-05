import pickle, json, sys
import heal_harness as H
from collections import defaultdict
sys.setrecursionlimit(2000000)
p=H.p
lin_rels=pickle.load(open('lin_rels.pkl','rb'))
N=H.NVARS; CONST=N  # value(CONST)=1
par=list(range(N+1)); mult=[1]*(N+1); off=[0]*(N+1)
def find2(x):
    if par[x]==x: return x,1,0
    r,m,o=find2(par[x])
    M=(mult[x]*m)%p; O=(mult[x]*o+off[x])%p
    par[x]=r; mult[x]=M; off[x]=O
    return r,M,O
contra=[]
def pin(v,k):  # x_v = k
    rv,mv,ov=find2(v)
    if rv==CONST:
        if (mv+ov-k)%p!=0: contra.append(('pin',v,k))
        return
    par[rv]=CONST; mult[rv]=0; off[rv]=((k-ov)*pow(mv,-1,p))%p
def union(a,b,m,o):  # x_a = m*x_b + o, m!=0
    ra,ma,oa=find2(a); rb,mb,ob=find2(b)
    if ra==rb:
        if (ma-(m*mb))%p!=0 or (oa-(m*ob+o))%p!=0: contra.append(('u',a,b))
        return
    if ra==CONST:  # x_a=ma+oa is constant K; then x_b=(K-o)/m
        K=(ma+oa)%p; pin(b,((K-o)*pow(m,-1,p))%p); return
    if rb==CONST:  # x_b=mb+ob=const L; x_a=m*L+o
        L=(mb+ob)%p; pin(a,(m*L+o)%p); return
    inv=pow((m*mb)%p,-1,p)
    par[rb]=ra; mult[rb]=(ma*inv)%p; off[rb]=((oa-m*ob-o)*inv)%p
multi=[]
for ai,rel,const in lin_rels:
    vs=list(rel)
    if len(vs)==1:
        v=vs[0]; c=rel[v]; pin(v,(-const*pow(c,-1,p))%p)
    elif len(vs)==2:
        va,vb=vs; ca,cb=rel[va],rel[vb]
        union(va,vb,(-cb*pow(ca,-1,p))%p,(-const*pow(ca,-1,p))%p)
    else:
        multi.append((rel,const))
print(f"unions done. contradictions: {len(contra)}; multi-term: {len(multi)}")
# reduce multi-term to roots
red=[]
for rel,const in multi:
    r=defaultdict(int); cc=const%p
    for v,c in rel.items():
        rv,mv,ov=find2(v)
        if rv==CONST: cc=(cc+c*(mv+ov))%p
        else: r[rv]=(r[rv]+c*mv)%p; cc=(cc+c*ov)%p
    r={k:v%p for k,v in r.items() if v%p}
    if r: red.append((r,cc))
    elif cc%p!=0: contra.append(('multi',))
print(f"reduced multi-term rels: {len(red)}; contradictions: {len(contra)}")
pickle.dump({'par':par,'mult':mult,'off':off,'red':red}, open('lin_reduced.pkl','wb'))
for a,b in [(14853,12186),(24908,16742)]:
    ra,ma,oa=find2(a); rb,mb,ob=find2(b)
    print(f"x_{a}:root{ra}(m{ma}) x_{b}:root{rb}(m{mb}) same={ra==rb} {'CONST-pinned' if ra==CONST or rb==CONST else ''}")
roots=set(r for rel,c in red for r in rel)
print(f"reduced system: {len(red)} eqs over {len(roots)} root-vars")
