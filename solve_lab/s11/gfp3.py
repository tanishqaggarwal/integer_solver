import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2
P=L.P
TARGETS=[('x3719',lambda v:v[3719]),('x25118',lambda v:v[25118]),
         ('x25614',lambda v:v[25614]),('x34220',lambda v:v[34220]),
         ('n-gap',lambda v:v[12186]-v[1308]),('m-gap',lambda v:v[24908]-v[19083])]
CTRL=[14515,19750,5096,21589,33708,31339,29261,26489,8060,19450,3473,8971,5616,245,
      27156,2467,19275,28548,6250,28486,5460,8363]
CFG={ '4bit':(542,47,438,91),
      '6bit-ab':(542,1685,47,1502,438,91),
      '6bit-cd':(542,47,490,438,1203,91),
      '8bit':(542,1685,47,1502,490,438,1203,91) }
for name,BITS in CFG.items():
    theta={c:0 for c in CTRL}
    v=engine2.close(BITS, theta, derive=False)
    r=[f(v)%P for _,f in TARGETS]
    J=[[0]*len(CTRL) for _ in TARGETS]
    for j,c in enumerate(CTRL):
        th=dict(theta); th[c]=1
        v1=engine2.close(BITS, th, derive=False)
        for i in range(len(TARGETS)): J[i][j]=(f_(v1)-f_(v))%P if (f_:=TARGETS[i][1]) else 0
    bad=fw.bad_checks(v)
    print(f"--- {name}: bad={len(bad)}")
    for i,(nm,_) in enumerate(TARGETS):
        print(f"    {nm:8s}: {[CTRL[j] for j in range(len(CTRL)) if J[i][j]]}")
    json.dump({'J':[[str(x) for x in row] for row in J],'r':[str(x) for x in r]}, open(f'J_{name}.json','w'))
