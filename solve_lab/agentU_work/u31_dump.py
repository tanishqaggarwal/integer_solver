import sys, json
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work')
import u20_sweep as S, umodel as U, uscore as SC, checker
v0=checker.load_assignment('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
for tag,(beta,a,b,src) in [('root_24601_2081',(U.ROOT,24601,2081,24601)),
                           ('root_generic_47',(U.ROOT,47,2081,47)),
                           ('slot28505_d6',(28505,438,2081,438))]:
    n,vv,sd=S.price(beta,a,b,src)
    d=sum(1 for i in range(38748) if vv[i]!=v0[i])
    p='u_build_%s.json'%tag
    json.dump({("x_%d"%i):vv[i] for i in range(38748) if vv[i]!=0}, open(p,'w'))
    print('%-18s -> %2d failing, %5d vars differ from deliverable -> %s'%(tag,n,d,p))
