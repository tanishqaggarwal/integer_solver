import pickle, sys, time
import heal_harness as H
from collections import defaultdict
p=H.p
R=pickle.load(open('lin_reduced.pkl','rb'))
par,mult,off,red=R['par'],R['mult'],R['off'],R['red']
CONST=H.NVARS
def find2(x):
    r=x; chain=[]
    while par[r]!=r: chain.append(r); r=par[r]
    m,o=1,0
    for v in reversed(chain):
        m=(mult[v]*m)%p; o=(mult[v]*o+off[v])%p  # wrong order; recompute properly
    # proper: value(x)=mult[x]*value(par[x])+off[x]; compose from root down
    # do it directly:
    M,O=1,0
    for v in chain:  # from x up to (not incl) root
        pass
    # simplest correct: recursive-free compose top-down
    M,O=1,0
    stack=chain
    # value(x) in terms of root: iterate x, par[x], ...: v_x = mult[x]*v_parx+off[x]
    # compose: start acc=(1,0) meaning identity on root; apply from root's child down to x
    acc_m,acc_o=1,0
    for v in reversed(chain):
        acc_m=(mult[v]*acc_m)%p; acc_o=(mult[v]*acc_o+off[v])%p
    return r,acc_m,acc_o
# core conditions: x_14853 ≡ x_12186, x_24908 ≡ x_16742
rows=[dict(rel) | {} for rel,c in red]  # copy
consts=[c for rel,c in red]
def rowfromvars(pairs, k):  # sum coef*x_v = k  -> reduce to roots
    r=defaultdict(int); cc=k
    for v,coef in pairs:
        rv,m,o=find2(v)
        if rv==CONST: cc=(cc-coef*(m+o))%p
        else: r[rv]=(r[rv]+coef*m)%p; cc=(cc - 0)%p; cc=(cc-coef*o)%p
    return {k2:v%p for k2,v in r.items() if v%p}, cc%p
for (a,b) in [(14853,12186),(24908,16742)]:
    r,cc=rowfromvars([(a,1),(b,-1)],0)  # x_a - x_b = 0
    if r: rows.append(r); consts.append(cc)
    elif cc%p!=0: print("CORE DIRECTLY INCONSISTENT",a,b)
print(f"total rows: {len(rows)}")
# sparse gaussian elimination mod p
col_rows=defaultdict(set)
for i,row in enumerate(rows):
    for v in row: col_rows[v].add(i)
alive=set(range(len(rows)))
pivots={}  # col -> row index
t0=time.time(); elim=0
# process by picking, among alive rows, pivot on the column with fewest alive rows
import heapq
def rowlen(i): return len(rows[i])
order=sorted(alive, key=rowlen)
contra=0
for i in list(order):
    if i not in alive: continue
    row=rows[i]
    if not row:
        if consts[i]%p!=0: contra+=1
        alive.discard(i); continue
    # pick pivot col = the var in row with fewest alive rows (Markowitz-ish)
    pc=min(row, key=lambda v: len(col_rows[v]&alive))
    pivots[pc]=i; alive.discard(i)
    inv=pow(row[pc],-1,p)
    # normalize
    rows[i]={v:(coef*inv)%p for v,coef in row.items()}; consts[i]=(consts[i]*inv)%p
    row=rows[i]
    # eliminate pc from all other alive rows containing it
    for j in list(col_rows[pc]):
        if j==i or j not in alive: continue
        f=rows[j].get(pc,0)
        if not f: continue
        for v,coef in row.items():
            nv=(rows[j].get(v,0)-f*coef)%p
            if nv: rows[j][v]=nv; col_rows[v].add(j)
            elif v in rows[j]: del rows[j][v]; col_rows[v].discard(j)
        consts[j]=(consts[j]-f*consts[i])%p
        if not rows[j] and consts[j]%p!=0: contra+=1; alive.discard(j)
    elim+=1
    if elim%1000==0: print(f"  eliminated {elim}, {time.time()-t0:.0f}s, contra={contra}",flush=True)
print(f"pivots: {len(pivots)}, contradictions: {contra}, time {time.time()-t0:.0f}s")
pickle.dump({'rows':rows,'consts':consts,'pivots':pivots}, open('gauss_result.pkl','wb'))
print("CONSISTENT" if contra==0 else f"INCONSISTENT ({contra} contradictions)")
