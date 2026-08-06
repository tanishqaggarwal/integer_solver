import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1
H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
H.val[9413]=0; H.val[17325]=0
H.forward()
H.val[8731]=H.val[4432]   # G2: x_19964 = x_4432
H.val[9118]=H.val[7068]   # G1: x_2099  = x_7068
H.forward()
F=H.fails()
G1=7376877*H.val[642]+H.val[2099]-H.val[7068]; G2=H.val[4432]-H.val[19964]-H.val[28730]
print('regime(1,1) clean: G1=%d G2=%d x_15298=%d fails=%d'%(G1,G2,H.val[15298],len(F)))
out={f'x_{i}':H.val[i] for i in range(H.NVARS) if H.val[i]!=0}
json.dump(out,open('sy_regime11_39018.json','w'))
print('saved sy_regime11_39018.json  nonzero vars=',len(out))
