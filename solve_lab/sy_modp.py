import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
p=H.p
B.regime11()
base=H.val[:]
F0=sorted(H.fails())
def resid_of(idxs):
    ns={'v':H.val,'__builtins__':{}}
    return {i:eval(H.eqcode[i],ns) for i in idxs}
R0=resid_of(F0)
def col(var,delta):
    H.val[var]=base[var]+delta; H.forward()
    R=resid_of(F0)
    H.val[var]=base[var]; H.forward()
    return {i:(R[i]-R0[i]) for i in F0}
# check affine: compare delta=1 vs delta=2 (scaled)
c1=col(8731,1); c1b=col(8731,2)
affine8731 = all(c1b[i]==2*c1[i] for i in F0)
c2=col(9118,1); c2b=col(9118,2)
affine9118 = all(c2b[i]==2*c2[i] for i in F0)
print('affine in x_8731?',affine8731,' in x_9118?',affine9118)
# build mod-p system: for each fail eq i: R0[i] + c1[i]*da + c2[i]*db ≡ 0 mod p
# solve 2-unknown system over GF(p) via least squares / consistency
rows=[(R0[i]%p, c1[i]%p, c2[i]%p) for i in F0]
# Solve: find da,db mod p s.t. r + a*da + b*db ≡0 for all rows. Use first two independent rows.
def inv(x): return pow(x%p,p-2,p)
# find 2 independent rows
import itertools
sol=None
for (r1,a1,b1),(r2,a2,b2) in itertools.combinations(rows,2):
    det=(a1*b2-a2*b1)%p
    if det!=0:
        da=((-r1)*b2 - (-r2)*b1)%p*inv(det)%p
        db=(a1*(-r2)-a2*(-r1))%p*inv(det)%p
        # verify all rows
        ok=all((r+a*da+b*db)%p==0 for r,a,b in rows)
        if ok:
            sol=(da,db); break
        else:
            sol=('inconsistent',(da,db)); 
print('mod-p solution (da,db) for x_8731,x_9118 shifts:', sol)
if sol and sol[0]!='inconsistent':
    da,db=sol
    print(' da=',da)
    print(' db=',db)
    # verify
    bad=[(r+a*da+b*db)%p for r,a,b in rows]
    print(' all rows ≡0 mod p?', all(x==0 for x in bad))
