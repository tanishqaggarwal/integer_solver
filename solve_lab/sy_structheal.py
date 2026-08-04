import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.val[4287]=1
H.val[31861]=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
H.val[14865]=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
H.val[8731]=r8; H.val[9118]=r9; H.val[9413]=0; H.val[17325]=0
H.forward()
H.val[4432]=H.val[19964]; H.val[7068]=H.val[2099]
H.forward()
H.val[950]=H.val[9106]//(13523997*p); H.val[6947]=(6122989*H.val[2239])//p; H.val[33168]=-(H.val[31731]//p)
H.forward()
print('before struct-heal:',len(H.fails()),'fails')
# heal residue-load atoms 7450, 7452 by setting free knobs EXACTLY
# 7450: x_2964 = x_26756 + x_579
# 7452: 9367949*(x_24548-x_25442) - p*x_11052 = 0 -> x_24548=x_25442, x_11052=0
H.val[2964]=H.val[26756]+H.val[579]
H.val[24548]=H.val[25442]; H.val[11052]=0
H.forward()
print('after healing 7450,7452:',len(H.fails()),'fails')
print(sorted(H.fails()))
