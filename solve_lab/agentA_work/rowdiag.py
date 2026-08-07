import sys, collections; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
from agrow import model
import amk_model as MK
A,K,R,rows,QUAD=model([37887,41906])
aff=[(e,c,lin) for e,c,lin,hq in rows if not hq]
nk=len(K); n=len(aff)
print('knobs',K)
print('%-8s %-6s %s'%('eq','nnz','knob support'))
for e,c,lin in aff:
    print('%-8d %-6d %s'%(e,len(lin),[K[j] for j in sorted(lin)]))
# column supports
col=collections.defaultdict(list)
for i,(e,c,lin) in enumerate(aff):
    for j in lin: col[j].append(e)
print()
for j,u in enumerate(K):
    print('col x%-6d support=%d %s'%(u,len(col[j]),col[j]))
