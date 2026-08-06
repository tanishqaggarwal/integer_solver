import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H
import sy_build as B
p=H.p
B.regime11()
def g(i): return H.val[i]
print('current loads:')
print(' x_9106 =',g(9106))
print(' x_2239 =',g(2239))
print(' x_31731=',g(31731))
print(' x_27177=',g(27177),' x_4306=',g(4306))
print('divisibility:')
print(' x_9106 % (13523997*p) ==0 ?', g(9106)%(13523997*p)==0, ' /p int?', g(9106)%p==0)
print(' x_2239 % p ==0 ?', g(2239)%p==0)
print(' x_31731 %p ==0?', g(31731)%p==0, ' ==0 exactly?', g(31731)==0)
# set handles optimally to kill p-parts
# 17897: x_9106 - 13523997*p*x_950 -> x_950 = x_9106//(13523997*p)
# 20866: 6122989*x_2239 - p*x_6947 -> x_6947 = 6122989*x_2239//p
# 20868: x_31731 + p*x_33168 -> x_33168 = -x_31731//p
H.val[950]=g(9106)//(13523997*p)
H.val[6947]=(6122989*g(2239))//p
H.val[33168]=-(g(31731)//p)
H.forward()
F=H.fails()
print('after optimal handles: fails=',len(F), sorted(F))
# residuals of remaining
ns={'v':H.val,'__builtins__':{}}
for i in sorted(F):
    r=eval(H.eqcode[i],ns)
    print(f'  eq {i}: resid%p={r%p}  resid/p~={r//p if abs(r)>p else r}  zeroModP={r%p==0}')
