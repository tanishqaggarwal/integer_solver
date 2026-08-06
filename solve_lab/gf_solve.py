"""Sparse Gaussian elimination over GF(p). Check consistency of J.delta = b and solve."""
import pickle, time, sys
p=2**256-2**32-977
SCR='/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad'

def rref_solve(rows, track_cert=False):
    """rows: list of (dict col->coef mod p, rhs mod p). Returns (consistent, pivots, echelon, cert).
    Reduces to echelon. If inconsistent, cert = list of (orig_row_idx, mult) giving 0=nonzero."""
    # echelon: pivotcol -> (rowdict, rhs, combo) where combo tracks original rows if track_cert
    piv={}  # col -> index into ech
    ech=[]  # list of [rowdict, rhs, combo]
    inconsistent=None
    for ri,(row0,rhs0) in enumerate(rows):
        row=dict(row0); rhs=rhs0%p
        combo={ri:1} if track_cert else None
        # reduce
        # iterate: while row has a col that is a pivot
        # process pivot cols in some order
        changed=True
        while True:
            pc=None
            for c in row:
                if c in piv:
                    pc=c; break
            if pc is None: break
            e=ech[piv[pc]]
            factor=row[pc]*e[0][pc]  # e[0][pc] is inverse-normalized? we normalize pivot to 1
            # e pivot normalized to 1 at pc, so factor = row[pc]
            factor=row[pc] % p
            # row -= factor * e
            er=e[0]
            for c,cf in er.items():
                nv=(row.get(c,0)-factor*cf)%p
                if nv: row[c]=nv
                elif c in row: del row[c]
            rhs=(rhs-factor*e[1])%p
            if track_cert:
                for k,mv in e[2].items():
                    combo[k]=(combo.get(k,0)-factor*mv)%p
        # now row has no pivot cols
        if not row:
            if rhs%p!=0:
                inconsistent=(ri, combo, rhs)
                # keep going? no, one certificate is enough
                return False, piv, ech, inconsistent
            # else dependent, skip
            continue
        # pick pivot col = smallest col (or first)
        pcol=min(row)
        inv=pow(row[pcol],p-2,p)
        nrow={c:(cf*inv)%p for c,cf in row.items()}
        nrhs=(rhs*inv)%p
        if track_cert:
            ncombo={k:(mv*inv)%p for k,mv in combo.items()}
        else:
            ncombo=None
        piv[pcol]=len(ech)
        ech.append([nrow,nrhs,ncombo])
    return True, piv, ech, None

def back_solve(piv, ech, ncols_free_default=0):
    """Return a particular solution dict col->value (free non-pivot cols = 0)."""
    x={}
    # process echelon rows in reverse (they're in creation order; pivot cols distinct)
    for pcol,idx in piv.items():
        pass
    # since non-pivot cols set to 0, x[pcol] = rhs - sum_{c!=pcol} row[c]*x[c]
    # need to solve in reverse dependency; echelon not sorted. Do iterative substitution.
    # Because each ech row is normalized (pcol coef=1), and other cols may be pivots of later rows.
    # Order by pivot col descending won't guarantee. Use: repeatedly evaluate rows whose non-pivot
    # unknown cols are all resolved or free(0). Simplest: set all non-pivot cols=0, then solve the
    # triangular-ish system by iterating to fixpoint.
    for _ in range(len(ech)+2):
        changed=False
        for pcol,idx in piv.items():
            row,rhs,_=ech[idx]
            s=rhs
            for c,cf in row.items():
                if c==pcol: continue
                s=(s-cf*x.get(c,0))%p
            if x.get(pcol,None)!=s%p:
                x[pcol]=s%p; changed=True
        if not changed: break
    return x

if __name__=='__main__':
    data=pickle.load(open(SCR+'/jac.pkl','rb'))
    Jrows=data['Jrows']; r=data['r']; breakable=data['breakable']
    # build row list: wiring rows (rhs 0) first, then residual rows (rhs -r)
    resid_atoms=[a for a in breakable if r.get(a,0)!=0]
    wiring=[a for a in breakable if a in Jrows and Jrows[a] and r.get(a,0)==0]
    rows=[]
    rowmeta=[]
    for a in wiring:
        rows.append((Jrows[a], 0)); rowmeta.append(a)
    for a in resid_atoms:
        rows.append((Jrows.get(a,{}), (-r[a])%p)); rowmeta.append(a)
    print(f"rows={len(rows)} (wiring={len(wiring)}, resid={len(resid_atoms)})")
    t0=time.time()
    ok,piv,ech,inc=rref_solve(rows, track_cert=True)
    print(f"GE done in {time.time()-t0:.1f}s. consistent={ok}. rank={len(piv)}")
    if not ok:
        ri,combo,rhs=inc
        print(f"INCONSISTENT at row {ri} (atom {rowmeta[ri]}), residual rhs={rhs}")
        # certificate: combo of original rows
        cert=[(rowmeta[k],mv) for k,mv in combo.items() if mv%p!=0]
        print(f"certificate size={len(cert)}")
        pickle.dump({'inc':inc,'combo':combo,'rowmeta':rowmeta,'cert':cert},open(SCR+'/cert.pkl','wb'))
    else:
        x=back_solve(piv,ech)
        # verify Jx=b
        bad=0
        for (row,rhs) in rows:
            s=sum(cf*x.get(c,0) for c,cf in row.items())%p
            if s!=rhs%p: bad+=1
        print(f"solution found; rows not matched={bad}")
        pickle.dump({'delta':x,'rowmeta':rowmeta},open(SCR+'/delta.pkl','wb'))
        print("saved delta.pkl")
