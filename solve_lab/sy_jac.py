import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
p=H.p
B.regime11()
base=H.val[:]                      # snapshot
def residuals():
    ns={'v':H.val,'__builtins__':{}}
    return [eval(c,ns) for c in H.eqcode]
R0=residuals()
F0=[i for i,r in enumerate(R0) if r!=0]
print('baseline fails:',len(F0), sorted(F0))
# composite moves: keep G1,G2 satisfied
# G2: x_4432 - x_19964(=x_8731) - p*x_9413 = 0
# G1: 7376877*p*x_17325 + x_2099(=x_9118) - x_7068 = 0
knobmoves={
 'h950':  {950:1},
 'h6947': {6947:1},
 'h33168':{33168:1},
 # move x_8731 & x_4432 together (keeps G2): affects x_19964 and loads and x_4432 roles
 'a8731': {8731:1, 4432:1},
 # move x_9118 & x_7068 together (keeps G1): affects x_2099 and loads and x_7068 roles
 'b9118': {9118:1, 7068:1},
 # p-granular G2 slack: x_9413 (x_28730=p*x_9413), and x_4432 += p to keep G2
 'g9413': {9413:1, 4432:p},
 # p-granular G1 slack: x_17325 (x_642=p*x_17325 -> 7376877*p), x_7068 += 7376877*p
 'g17325':{17325:1, 7068:7376877*p},
}
def effect(mv):
    for k,v in mv.items(): H.val[k]=base[k]+v
    H.forward()
    R=residuals()
    d={i:(R[i]-R0[i]) for i in range(len(R)) if R[i]!=R0[i]}
    for k,v in mv.items(): H.val[k]=base[k]
    H.forward()
    return d
cols={}
for name,mv in knobmoves.items():
    d=effect(mv)
    cols[name]=d
    aff=sorted(d.keys())
    outside=[e for e in aff if e not in F0]
    print(f'{name}: affects {len(aff)} eqs, {len(outside)} OUTSIDE fails: {outside[:12]}')
# Save R0 and cols
import pickle
pickle.dump({'R0':R0,'F0':F0,'cols':cols,'knobmoves':knobmoves}, open('sy_jac.pkl','wb'))
