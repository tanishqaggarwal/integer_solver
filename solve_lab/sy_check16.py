import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
import sy_build as B
p=H.p
r8=109044024666698959972204451600908701898659086097062528124234304603594878834481
r9=33371159155735472537534252650716501592825364489306217536352743247010353604716
# agentA free-leaf values
vA=H.loadd('best_agentA_39022.json')
V4432=vA.get(4432,0); V7068=vA.get(7068,0)
print('agentA x_4432 %p =', V4432%p)
print('r8            =', r8)
print('match r8? ', V4432%p==r8)
print('agentA x_7068 %p =', V7068%p, ' r9=',r9,' match?', V7068%p==r9)
# examine the 16 outside eqs: load atoms structure. Check what breaks when x_4432 -> large.
# Build regime11, then move ONLY x_4432 by p (via x_9413) and see which outside eqs break + their atoms
B.regime11()
base=H.val[:]
F0=set(H.fails())
# move x_4432 by +p (x_9413 +=1), x_7068 by +7376877*p (x_17325+=1) to keep G1/G2, see new fails
H.val[9413]=1; H.val[4432]=base[4432]+p
H.forward()
newF=set(H.fails())
broke=sorted(newF-F0)
print('moving x_4432 by p (x_9413=1): new broken eqs=',broke)
