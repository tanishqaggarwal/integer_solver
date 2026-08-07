"""W4: exact scoring harness.  Diagnose the deliverable's lie, then try to repair it
with a PIN lie instead of a ROUTE lie.  Every number here comes from checker.py's own
compiled equations -- no incidence anywhere."""
import sys, json, pickle, collections, math, time
sys.path.insert(0,'/home/user/integer_solver/solve_lab')
import checker
B='/home/user/integer_solver/solve_lab/agentU_work/'
L=pickle.load(open(B+'v_leaves.pkl','rb')); XY=pickle.load(open(B+'w_xy.pkl','rb'))
D=pickle.load(open(B+'v_defs.pkl','rb'))
X=XY['X']; Y=XY['Y']; sel2exp=L['sel2exp']; p=L['p']
t0=time.time(); codes,varsets=checker.load_equations(); print('equations loaded %.1fs'%(time.time()-t0))
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
def score(v):
    f=checker.evaluate_all(codes,v); return len(f),f
n,f=score(v0); print('deliverable: %d failing %s'%(n,f))
ea,eb=72,235
print('leaf 72 : Xwire x%d M=%d z=x%d ; Ywire x%d M=%d z=x%d'%(X[ea][0],X[ea][1],X[ea][2],Y[ea][0],Y[ea][1],Y[ea][2]))
print('leaf 235: Xwire x%d M=%d z=x%d ; Ywire x%d M=%d z=x%d'%(X[eb][0],X[eb][1],X[eb][2],Y[eb][0],Y[eb][1],Y[eb][2]))
for e in (ea,eb):
    print(' e=%d  deliverable X wire=%s  rawC=%s'%(e,str(v0[X[e][0]])[:24],str(X[e][3])[:24]))
    print(' e=%d  deliverable Y wire=%s  rawC=%s'%(e,str(v0[Y[e][0]])[:24],str(Y[e][3])[:24]))
# what values does the deliverable put on wires?  how many carry leaf72 vs leaf235 coords
vals=collections.Counter(v0[i] for i in range(len(v0)) if v0[i]!=0)
print('distinct non-zero values in the deliverable:',len(vals),' most common:',[(str(k)[:14],c) for k,c in vals.most_common(4)])
tag={}
for e in (ea,eb):
    tag[('X',e)]=X[e][3]%p; tag[('Y',e)]=Y[e][3]%p
for k,c in vals.most_common(8):
    m=[kk for kk,vv in tag.items() if vv==c%p]
    if m: print('   value carried by %d wires == %s of leaf %d'%(vals[k],m[0][0],m[0][1]))
# --- CANDIDATE: replace the ROUTE lie with a PIN lie on leaf 235 ---
res=[]
for src,dst in ((ea,eb),(eb,ea)):
    v=list(v0)
    okX=(X[src][3]-X[dst][3])%X[dst][1]==0
    zx=(X[src][3]-X[dst][3])//X[dst][1] if okX else None
    v[X[dst][0]]=X[src][3]; v[Y[dst][0]]=Y[src][3]
    if okX: v[X[dst][2]]=zx
    v[Y[dst][2]]=Y[src][3]-Y[dst][3]
    n2,f2=score(v)
    res.append(('pin lie: leaf %d wires := leaf %d constants'%(dst,src), okX, n2))
    print('  %-46s X-divisible=%-5s -> %d failing'%(res[-1][0],okX,n2))
# --- CANDIDATE: joint CRT value W on both leaves ---
Ma,Mb=X[ea][1],X[eb][1]; g=math.gcd(Ma,Mb); diff=X[ea][3]-X[eb][3]
print('joint CRT: gcd(Ma,Mb)=%d  divides diff: %s'%(g,diff%g==0))
if diff%g==0:
    lcm=Ma//g*Mb
    # W = C_aX + Ma*t  with  C_aX + Ma*t = C_bX (mod Mb)
    t=((X[eb][3]-X[ea][3])//g * pow((Ma//g)%(Mb//g), -1, Mb//g)) % (Mb//g)
    W=X[ea][3]+Ma*t
    assert (W-X[ea][3])%Ma==0 and (W-X[eb][3])%Mb==0
    WY=Y[ea][3]
    v=list(v0)
    for e in (ea,eb):
        v[X[e][0]]=W; v[X[e][2]]=(W-X[e][3])//X[e][1]
        v[Y[e][0]]=WY; v[Y[e][2]]=WY-Y[e][3]
    n3,f3=score(v); print('  joint CRT common W on BOTH leaves           -> %d failing'%n3)
    v2=list(v); n4,_=score(v2)
pickle.dump({'codes_n':len(codes)}, open(B+'w_score_meta.pkl','wb'))
