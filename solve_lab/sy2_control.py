import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
fc=H.loadd('fc_partial.json')
def setup():
    for v in H.freeinp: H.val[v]=fc.get(v,0)
    H.forward()
setup()
F0=set(H.fails())
print('fc_partial fails:',len(F0))
# pin handles: check free
# 3277: x_2081*(x_6418-K1)-15804267*x_26777  -> x_26777 handle for x_6418
# 3279: x_2081*(x_12553-K2)-x_13458          -> x_13458 for x_12553
# 3269: x_4287*(x_31861-K3)-13479571*x_27676 -> x_27676 for x_31861
# 3271: x_4287*(x_14865-K4)-x_7574           -> x_7574 for x_14865
handles={6418:26777,12553:13458,31861:27676,14865:7574}
for pin,hnd in handles.items():
    print('pin x_%d (%s) handle x_%d (%s)'%(pin,'FREE' if pin in H.freeinp else 'gate',hnd,'FREE' if hnd in H.freeinp else 'gate'))
# candidate knobs to control x_27177,x_4306 mod p: the pins directly (they're free), and x_8731,x_9118
cands=[6418,12553,31861,14865,8731,9118]
def state():
    return H.val[27177]%p, H.val[4306]%p
setup(); u0,w0=state()
print('baseline x_27177%p=',u0,' x_4306%p=',w0)
for kn in cands:
    setup()
    old=H.val[kn]; H.val[kn]=old+1; H.forward()
    u1,w1=state()
    # ripple: broken eqs outside F0
    nb=set(H.fails())-F0
    du=(u1-u0)%p; dw=(w1-w0)%p
    setup()
    print('knob x_%d: du27177=%s dw4306=%s  ripple_new=%d'%(kn,'0' if du==0 else 'nz',('0' if dw==0 else 'nz'),len(nb)))
