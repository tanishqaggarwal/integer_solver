import pickle, re, collections, sys
m=pickle.load(open('model3.pkl','rb')); atoms=m['atoms']
d=pickle.load(open('dag.pkl','rb')); info=d['info']; free=d['free']
dd=pickle.load(open('definer.pkl','rb')); definer=dd['definer']
p2=pickle.load(open('prop2.pkl','rb')); v=p2['v']
avars=[sorted(vs) for _,vs in info]
aid=int(sys.argv[1])
seen=set(); stack=list(avars[aid]); order=[]
while stack:
    u=stack.pop()
    if u in seen: continue
    seen.add(u); order.append(u)
    di=definer[u]
    if di is not None:
        for w in avars[di]:
            if w!=u: stack.append(w)
print("ATOM",aid,atoms[aid])
for u in sorted(seen):
    di=definer[u]
    print(f"  x_{u} = {str(v[u])[:70]:70s}  def: {atoms[di][:90] if di is not None else 'FREE'}")
