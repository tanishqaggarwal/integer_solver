"""Refined construction: activate branch, fix the newly-lit pins deliberately, measure."""
import pickle, sys
import harness as H
exec(open('repair.py').read().split('if __name__')[0])
P=2**256-2**32-977
roots=pickle.load(open('roots.pkl','rb'))
checks=[a for a in range(len(polys)) if a not in atom_out]
rp={a:(roots[a] if a in roots else polys[a]) for a in checks}
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319
H1=111427306269069086410860162925174701533202290689960545658556486068216322467769731977368545
H2=80924611027189446969163401728848115572253569416947377083869986344878002353122806121600425

def nz(v): return sorted(a for a,Pp in rp.items() if evalpoly(Pp,v)!=0)

def build(b=2527):
    v=H.load_assignment('../best/new_instance_partial_39022.json')
    ripple(v,{b:1, 5096:K1, 33612:0})
    ripple(v,{14853:v[12186]})
    ripple(v,{7068:v[2099]+7376877*v[642], 4432:v[19964]+v[28730]})
    ripple(v,{24548:v[25442]})
    print('A. after branch activation      :', nz(v))
    ripple(v,{20742:H1, 8824:H2})           # satisfy the pins the bit lit up
    print('B. after satisfying lit pins    :', nz(v))
    ripple(v,{16742:v[19083]})              # restore x_16742 pin (atom 26731)
    ripple(v,{7068:v[2099]+7376877*v[642], 4432:v[19964]+v[28730]})
    ripple(v,{24548:v[25442], 14853:v[1308]})
    print('C. after re-aligning mirrors    :', nz(v))
    return v

if __name__=='__main__':
    b=int(sys.argv[1]) if len(sys.argv)>1 else 2527
    v=build(b)
    codes,_=H.load_equations(); f=H.evaluate(codes,v)
    print(f'EQUATIONS satisfied: {len(codes)-len(f)}/{len(codes)}  ({len(f)} failing)')
    H.save_assignment(v, f'construct2_{b}.json')
