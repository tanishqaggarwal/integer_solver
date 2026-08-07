import sys, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentU_work/mirror')
import engine3 as E3, harness as H
eng=E3.Eng(E3.BASE_DEMOTE)
S=pickle.load(open('x_seed.pkl','rb')); sd=S['seed']
XY=pickle.load(open('w_xy.pkl','rb')); SU=pickle.load(open('v_supp2.pkl','rb'))
L=pickle.load(open('v_leaves.pkl','rb'))
supp=SU['supp']
CX=XY['X'][72][3]; CY=XY['Y'][72][3]
sel2exp=L['sel2exp']
def sname(u):
    s=supp.get(u)
    if s is None: return 'no-supp'
    ex=sorted(sel2exp[t] for t in s if t in sel2exp)
    return '|supp|=%d %s'%(len(ex), ex if len(ex)<=6 else str(ex[:3])+'...'+str(ex[-3:]))
print('== BASE_DEMOTE atoms and the vars they define ==')
for a in E3.BASE_DEMOTE:
    print('  atom %d: %s   -> defines x_%s'%(a,H.atoms[a],E3.ATOM2VAR.get(a)))
print()
print('== route wires (seed value == leaf72 X or Y const) ==')
rx=[k for k in sd if sd[k]==CX]; ry=[k for k in sd if sd[k]==CY]
print(' X-route (%d): %s'%(len(rx),sorted(rx)))
for k in sorted(rx): print('    x_%-6d %s'%(k,sname(k)))
print(' Y-route (%d): %s'%(len(ry),sorted(ry)))
for k in sorted(ry): print('    x_%-6d %s'%(k,sname(k)))
print()
print('== remaining seed entries ==')
for k in sorted(sd):
    if k in rx or k in ry: continue
    print('  x_%-6d bits=%-5d %s'%(k,sd[k].bit_length(),sname(k)))
