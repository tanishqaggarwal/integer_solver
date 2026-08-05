"""Compact GF(p) Newton heal: from a config, iteratively move free inputs (least-norm mod p)
to zero the currently-failing eqs. Confirms whether continuous healing can beat 11."""
import _bitlab as L, heal_harness as H
from jac_lib import D, freeidx, freelist
p=H.p
def inv(a): return pow(a%p,p-2,p)

def heal(startfn, iters=8, tag=''):
    startfn()
    for it in range(iters):
        F=H.fails()
        if len(F)==0:
            print(f'  [{tag}] iter{it}: SOLVED 0 fails!'); return 0
        # exact GF(p) jacobian of failing eqs wrt all free inputs
        vd=[None]*H.NVARS
        for j in H.freeinp: vd[j]=D(H.val[j],{freeidx[j]:1})
        ns={'v':vd,'__builtins__':{}}
        for k,t in enumerate(H.order):
            r=eval(H.gcode[k],ns); vd[t]=r if isinstance(r,D) else D(r)
        rows=[]; rhs=[]
        for i in F:
            rr=eval(H.eqcode[i],{'v':vd,'__builtins__':{}})
            if isinstance(rr,D): rows.append(dict(rr.g)); rhs.append((-rr.v)%p)
            else: rows.append({}); rhs.append((-(rr%p))%p)
        # least-norm-ish solve J d = rhs mod p via row reduction; then apply move to free inputs
        cols=sorted(set(c for g in rows for c in g))
        # gaussian elim to get a particular solution d (col-> value)
        R=[dict(g) for g in rows]; B=list(rhs); nr=len(R)
        used=[False]*nr; pivot={}
        for c in cols:
            pr=-1
            for r in range(nr):
                if not used[r] and R[r].get(c,0)!=0: pr=r; break
            if pr<0: continue
            used[pr]=True; iv=inv(R[pr][c])
            R[pr]={k:(v*iv)%p for k,v in R[pr].items()}; B[pr]=(B[pr]*iv)%p
            for r in range(nr):
                if r!=pr and R[r].get(c,0)!=0:
                    f=R[r][c]
                    for k,v in R[pr].items():
                        nv=(R[r].get(k,0)-f*v)%p
                        if nv: R[r][k]=nv
                        elif k in R[r]: del R[r][k]
                    B[r]=(B[r]-f*B[pr])%p
            pivot[c]=pr
        incon=sum(1 for r in range(nr) if not R[r] and B[r]!=0)
        # apply particular solution: free var col c -> move by B[pivot[c]]
        d={}
        for c,pr in pivot.items(): d[c]=B[pr]
        for c,val in d.items():
            fv=freelist[c]
            H.val[fv]=(H.val[fv]+val)%p  # move mod p
        H.forward()
        print(f'  [{tag}] iter{it}: fails={len(F)}, jac-inconsistent-rows={incon}, moved {len(d)} free inputs')
    return len(H.fails())

def start_close():
    L.apply_pattern(L.AGENTA_BITS,twopass=False)
    V=H.val; V[7068]=V[2099]+7376877*V[642]; V[4432]=V[19964]+V[28730]; H.forward()

print('Newton heal from close-both (start 16):')
r=heal(start_close, iters=6, tag='close16')
print('final fails:', r)
