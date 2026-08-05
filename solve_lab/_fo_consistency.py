import pickle, time
from collections import defaultdict
p=2**256-2**32-977
D=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac.pkl','rb'))
rows=[dict(r) for r in D['rows']]; consts=list(D['consts'])
n=len(rows)
col_rows=defaultdict(set)
for i,row in enumerate(rows):
    for v in row: col_rows[v].add(i)
alive=set(range(n))
pivots={}; contra=0; t0=time.time(); elim=0
order=sorted(range(n), key=lambda i: len(rows[i]))
for i in order:
    if i not in alive: continue
    row=rows[i]
    if not row:
        if consts[i]%p!=0: contra+=1
        alive.discard(i); continue
    pc=min(row, key=lambda v: len(col_rows[v]&alive))
    pivots[pc]=i; alive.discard(i)
    inv=pow(row[pc],-1,p)
    rows[i]={v:(c*inv)%p for v,c in row.items()}; consts[i]=(consts[i]*inv)%p
    row=rows[i]
    for j in list(col_rows[pc]):
        if j==i or j not in alive: continue
        f=rows[j].get(pc,0)
        if not f: continue
        for v,c in row.items():
            nv=(rows[j].get(v,0)-f*c)%p
            if nv: rows[j][v]=nv; col_rows[v].add(j)
            elif v in rows[j]: del rows[j][v]; col_rows[v].discard(j)
        consts[j]=(consts[j]-f*consts[i])%p
        if not rows[j] and consts[j]%p!=0: contra+=1; alive.discard(j)
    elim+=1
print(f"pivots(rank J)={len(pivots)}  contradictions={contra}  time {time.time()-t0:.1f}s")
print("FIRST-ORDER NEWTON:", "CONSISTENT (delta exists mod p)" if contra==0 else f"INCONSISTENT mod p ({contra} contradiction rows)")
