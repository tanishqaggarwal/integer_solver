import os,sys,json
os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
p=H.p
vA=H.loadd('best_agentA_39022.json')
for v in H.freeinp: H.val[v]=vA.get(v,0)
H.forward()
F0=set(H.fails())
print("baseline fails:",sorted(F0))
# current gate values of x_2099, x_19964
x2099=H.val[2099]; x19964=H.val[19964]
print("x_2099=",x2099%p," x_19964=",x19964%p)
# Set x_7068 = x_2099, x_4432 = x_19964 (free inputs), keep x_17325=x_9413=0
H.val[7068]=x2099
H.val[4432]=x19964
H.forward()  # recompute gates
F1=set(H.fails())
print("\nAfter x_7068<-x_2099, x_4432<-x_19964:")
print("  #fails:",len(F1))
print("  fixed:",sorted(F0-F1))
print("  NEW broken:",sorted(F1-F0))
print("  still fail:",sorted(F1&F0))
