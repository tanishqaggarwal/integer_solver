"""U24: does price(beta,a,b,src=a) depend on a at all, or only on (beta, the lying leaf b)?"""
import sys, time, collections
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S, umodel as U
# pick a slot with several leaves on each side
cand=[(b,len(U.LIVELEAF[U.tree[b][0]]),len(U.LIVELEAF[U.tree[b][1]])) for b in U.SLOTS]
cand=[c for c in cand if 3<=c[1]<=6 and 3<=c[2]<=6]
for beta,ni,nj in cand[:3]:
    A=sorted(U.LIVELEAF[U.tree[beta][0]]); B=sorted(U.LIVELEAF[U.tree[beta][1]])
    print('== beta=%d |I|=%d |J|=%d'%(beta,ni,nj))
    tab=collections.defaultdict(set)
    for a in A:
        row=[]
        for b in B:
            n,_,_=S.price(beta,a,b,a)   # b's chain lies
            row.append(n); tab[b].add(n)
        print('   a=%-6d  src=a, over b: %s'%(a,row))
    print('   -> price depends only on lying leaf b? ', all(len(v)==1 for v in tab.values()))
