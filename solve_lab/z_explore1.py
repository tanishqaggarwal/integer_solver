import os,sys
os.chdir('/home/user/integer_solver/solve_lab')
sys.path.insert(0,'.')
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F=H.fails()
print("fails:",F)
kv=[7068,2099,642,17325,4432,19964,28730,9413,2081,4287,
    31861,6418,9118,14865,12553,8731,17499,28599,26064]
print("--- key var values ---")
for v in kv:
    val=H.val[v]; r=val%p
    tag=''
    if val==p: tag='==p'
    elif val==0: tag='==0'
    elif val==p*(val//p) and val!=0: tag=f'=={val//p}*p'
    print(f"x_{v}={val}"[:70],f"  mod p={r}  {tag}")
