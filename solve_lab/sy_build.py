import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
K3=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110
K4=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706
def base():
    vA=H.loadd('best_agentA_39022.json')
    for v in H.freeinp: H.val[v]=vA.get(v,0)
    return vA
def regime11(extra=None):
    base()
    H.val[4287]=1          # x_9062=1
    H.val[31861]=K3        # satisfy pin atom 3269
    H.val[14865]=K4        # satisfy pin atom 3271
    H.val[9413]=0; H.val[17325]=0
    H.forward()
    # absorb G1,G2 (free knobs x_8731, x_9118)
    H.val[8731]=H.val[4432]   # x_19964=x_4432 -> G2=0
    H.val[9118]=H.val[7068]   # x_2099 =x_7068 -> G1=0
    if extra:
        for k,v in extra.items(): H.val[k]=v
    H.forward()
if __name__=='__main__':
    regime11()
    F=H.fails()
    G2=H.val[4432]-H.val[19964]-H.val[28730]
    G1=7376877*H.val[642]+H.val[2099]-H.val[7068]
    print('G1=',G1,'G2=',G2,'x_15298=',H.val[15298])
    print('FAILS=',len(F))
    print(sorted(F))
