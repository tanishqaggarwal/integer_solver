import heal_harness as H, jac_lib as J, pickle
p=H.p
D=pickle.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac.pkl','rb'))
rows=[dict(r) for r in D['rows']]; consts=list(D['consts'])
F=[2554, 6816, 8124, 8680, 9421, 12231, 12270, 12350, 14584, 22044, 29125]; Fset=set(F)
def rank_of(rowset, withconst=False):
    piv={}   # col->rowdict ; if withconst, rowdict includes special key 'C'
    contra=0
    for i in rowset:
        rd=dict(rows[i])
        if withconst: rd['C']=consts[i]%p
        # reduce
        while True:
            cols=[c for c in rd if c!='C' and rd[c]%p]
            if not cols: break
            c=min(cols)
            if c in piv:
                prow=piv[c]; f=rd[c]
                for k,v in prow.items():
                    nv=(rd.get(k,0)-f*v)%p
                    if nv: rd[k]=nv
                    elif k in rd: del rd[k]
            else:
                inv=pow(rd[c],-1,p); piv[c]={k:(v*inv)%p for k,v in rd.items()}; break
        else:
            pass
        # after loop, check contradiction (only const left)
        realcols=[c for c in rd if c!='C' and rd[c]%p]
        if withconst and not realcols and rd.get('C',0)%p!=0: contra+=1
    return len(piv),contra
allrows=list(range(len(rows)))
sat=[i for i in allrows if i not in Fset]
rs,_=rank_of(sat)
ra,_=rank_of(allrows)
rc,contra=rank_of(allrows,withconst=True)
print(f"rank(J_sat)={rs}")
print(f"rank(J_all)={ra}")
print(f"rank([J_all|-r])={rc}  contradictions during reduce={contra}")
print(f"fails add {ra-rs} to rank; augmented adds {rc-ra} => {'INCONSISTENT' if rc>ra else 'CONSISTENT'}")
print(f"number of independent first-order obstructions = {rc-ra}")
