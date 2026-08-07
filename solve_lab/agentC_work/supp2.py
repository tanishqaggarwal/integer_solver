import sys, collections, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/s9/eff')
import lib as L
NV=L.NVARS
outs={}
for a,(c,t) in L.atom_out.items(): outs[t]=a
free=[v for v in range(NV) if v not in outs]
freeset=set(free)
# graph: v -> inputs
adj={}
for v in range(NV):
    if v in outs:
        adj[v]=[u for u in L.avars[outs[v]] if u!=v]
    else: adj[v]=[]
# iterative Tarjan SCC
index=[0]; idx={}; low={}; onstk={}; stk=[]; sccs=[]
for root in range(NV):
    if root in idx: continue
    work=[(root,0)]
    while work:
        v,pi=work[-1]
        if pi==0:
            idx[v]=low[v]=index[0]; index[0]+=1; stk.append(v); onstk[v]=True
        rec=False
        for i in range(pi,len(adj[v])):
            w=adj[v][i]
            if w not in idx:
                work[-1]=(v,i+1); work.append((w,0)); rec=True; break
            elif onstk.get(w): low[v]=min(low[v],idx[w])
        if rec: continue
        if low[v]==idx[v]:
            c=[]
            while True:
                w=stk.pop(); onstk[w]=False; c.append(w)
                if w==v: break
            sccs.append(c)
        work.pop()
        if work:
            u=work[-1][0]; low[u]=min(low[u],low[v])
sizes=collections.Counter(len(c) for c in sccs)
print('SCC size hist:',sorted(sizes.items())[:10], 'max', max(len(c) for c in sccs))
big=[c for c in sccs if len(c)>1]
print('nontrivial SCCs:',len(big),'total vars in them',sum(len(c) for c in big))
pickle.dump(sccs,open('/home/user/integer_solver/solve_lab/agentC_work/sccs.pkl','wb'))
for c in sorted(big,key=len,reverse=True)[:5]:
    print('SCC size',len(c),'sample',c[:6])
    for v in c[:3]:
        print('   x_%d <- a%d: %s'%(v,outs[v],L.atom_src[outs[v]][:130]))
