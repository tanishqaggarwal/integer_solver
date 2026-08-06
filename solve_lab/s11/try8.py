import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 's9', 'eff'))
import lib as L, fw, engine2, deep
B4=(542,47,438,91)
B8=(542,1685, 47,1502, 490,438, 1203,91)
for name,BITS in [("4bit",B4),("8bit",B8)]:
    t0=time.time()
    v=engine2.close(BITS, {})
    b=fw.bad_checks(v); av=L.all_atom_values(v); f=L.failing_eqs(av)
    print(f"{name}: bad={len(b)} failing={len(f)} score={L.NEQ-len(f)} ({time.time()-t0:.0f}s)")
    print("   bad:", b)
    if len(b)<=10:
        locked=set(BITS)|engine2.DERIVED
        for a in b:
            h,base=deep.handles(v,a,locked=locked)
            modp=[(t,d) for t,d in h if d%L.P!=0]
            print(f"   a{a}: handles={len(h)} MODP-live={[(t,len(L.var_atoms[t])) for t,d in modp]}")
