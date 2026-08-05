import heal_harness as H, jac_lib as J, time, pickle
p=H.p
d=H.loadd('best_agentA_39022.json')
base={v:d.get(v,0) for v in H.freeinp}
for v in H.freeinp: H.val[v]=base[v]
H.forward()
t0=time.time()
vd=J.build_duals()
print(f"duals built {time.time()-t0:.1f}s",flush=True)
NEQ=len(H.eqcode)
rows=[]      # list of dict col->coef (mod p)
consts=[]    # -r_i  (RHS for J*delta = -r)
t0=time.time()
nz_resid=0
for i in range(NEQ):
    rv,grad=J.eq_jac_row(i,vd)
    rows.append(grad)
    consts.append((-rv)%p)
    if rv!=0: nz_resid+=1
    if (i+1)%8000==0: print(f"  {i+1}/{NEQ} eqs, {time.time()-t0:.1f}s",flush=True)
print(f"jacobian rows built {time.time()-t0:.1f}s; eqs with nonzero resid mod p = {nz_resid}",flush=True)
# stats
tot=sum(len(r) for r in rows)
print(f"total nonzeros in J = {tot}; avg row len {tot/NEQ:.1f}")
# columns actually used
usedcols=set()
for r in rows: usedcols|=set(r.keys())
print(f"free-input columns with any nonzero Jacobian entry = {len(usedcols)} / {J.NF}")
pickle.dump({'rows':rows,'consts':consts,'nz_resid':nz_resid},open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/jac.pkl','wb'))
print("saved jac.pkl")
