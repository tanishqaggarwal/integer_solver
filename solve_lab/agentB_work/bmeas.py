import os,sys,json,collections
os.environ['ORIENT']=os.environ.get('ORIENT','orient7.pkl'); sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentB_work')
import beval as E, bfix as F
P=F.P; Q=F.Q
v0=E.load('/home/user/integer_solver/solve_lab/best/new_instance_partial_39026.json')
print('x15616=',v0[15616]==P and 'p' or v0[15616].bit_length())
print('x13859=',v0[13859]==P and 'p' or v0[13859].bit_length())
print('x22665=',v0[22665]==P and 'p' or v0[22665].bit_length())
def ev(shift):
    fv={v:v0[v] for v in E.free}
    for k,d in shift.items(): fv[k]=fv[k]+d
    val,nd,_=E.forward(fv,default=v0)
    return val
base=ev({})
def probe(knob, d):
    val=ev({knob:d})
    return dict(
      A=(val[28730]%P),
      B=(val[8731]%P),
      Ee=(val[9118]%P),
      D=((val[7068]-val[2099])%Q),
      G1=(val[7927]%v0[15616]),
      G2=((val[15324]-val[37254])%8481759),
      G3=(val[579]%v0[13859]),
    )
b=probe(4432,0)
print('base residues:',{k:(v.bit_length() if v else 0) for k,v in b.items()})
for knob,d,lab in [(4432,1,'x4432+1'),(4432,P,'x4432+p'),(7068,1,'x7068+1'),(7068,Q,'x7068+Q'),(3349,P,'x3349+p'),(9118,P,'x9118+p')]:
    r=probe(knob,d)
    print(lab, {k:('same' if r[k]==b[k] else 'MOVED') for k in r})
# linear derivative of G1,G2,G3 wrt x4432, x7068
for knob in [4432,7068]:
    for dd in [1]:
        r1=probe(knob,dd); r2=probe(knob,2*dd)
        print('knob x%d: dG1=%s dG2=%s dG3=%s'%(knob,
          (r1['G1']-b['G1']), (r1['G2']-b['G2']), (r1['G3']-b['G3'])))
