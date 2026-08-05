import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.forward()
print('fc_partial fails:',len(H.fails()))
def info(idx):
    tag='FREE' if idx in H.freeinp else 'gate'
    print('  x_%d: %s val=%s'%(idx,tag,H.val[idx]))
print('--- sinks: free or gate? ---')
for idx in [9629,30095,950,23754,26874,6947,37720,8976,33168,35619,24490]:
    info(idx)
print('--- gadget load values ---')
for idx in [9106,2239,31731,27177,4306]:
    print('  x_%d=%s'%(idx,H.val[idx]))
print('x_9106 mod 13523997 =',H.val[9106]%13523997)
print('x_9106 mod p =',H.val[9106]%p)
print('x_9106 mod (13523997*p) =',H.val[9106]%(13523997*p))
print('x_31731 ==0?',H.val[31731]==0,' mod p =',H.val[31731]%p)
