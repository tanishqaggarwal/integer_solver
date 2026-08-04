import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.forward()
print('base fails:',len(H.fails()))
# core move: x_14865=x_12553, x_31861=x_6418  => x_27019=x_17925=0 => loads vanish
H.val[14865]=H.val[12553]
H.val[31861]=H.val[6418]
H.val[33168]=0
H.forward()
print('x_27019=',H.val[27019],'x_17925=',H.val[17925])
print('x_27177=',H.val[27177],'x_4306=',H.val[4306])
print('x_31731=',H.val[31731],'x_9106=',H.val[9106],'x_2239=',H.val[2239])
G1=7376877*H.val[642]+H.val[2099]-H.val[7068]; G2=H.val[4432]-H.val[19964]-H.val[28730]
print('G1=',G1,'G2=',G2,'x_15298=',H.val[15298])
F=H.fails()
print('fails after core move:',len(F))
print(sorted(F))
