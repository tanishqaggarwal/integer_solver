"""Independent test of the elliptic-curve hypothesis, and of whether my two congruence
residues D0 = (x7068-x2099) mod p and K2 = x28730 mod p live on that curve."""
import sys, json; sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentA_work')
import env, lib as L
P=env.P
SEC_B=7
N_SEC=0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
COORD={'x1':12186,'y1':16742,'x2':14853,'y2':24908,'x3':22162,'y3':30213}
def qr(a):
    a%=P
    return a==0 or pow(a,(P-1)//2,P)==1
for path in sys.argv[1:]:
    v=L.load(path); av=L.all_atom_values(v)
    s=L.NEQ-len(L.failing_eqs(av))
    print('=== %s (score %d) ==='%(path.split('/')[-1],s))
    C={k:v[u]%P for k,u in COORD.items()}
    for k in ['x1','y1','x2','y2','x3','y3']:
        print('   %-3s = x%-6d  mod p = %d'%(k,COORD[k],C[k]))
    for (xk,yk) in [('x1','y1'),('x2','y2'),('x3','y3')]:
        b=(C[yk]*C[yk]-pow(C[xk],3,P))%P
        print('   b(%s,%s) = y^2 - x^3 = %d %s'%(xk,yk,b,'  <== 7 : secp256k1' if b==7 else ''))
    b1=(C['y1']**2-pow(C['x1'],3,P))%P
    b2=(C['y2']**2-pow(C['x2'],3,P))%P
    b3=(C['y3']**2-pow(C['x3'],3,P))%P
    print('   COMMON CURVE?  b1==b2: %s   b1==b3: %s'%(b1==b2,b1==b3))
    # affine addition check (x3,y3) = (x1,y1)+(x2,y2) on b1
    if C['x1']!=C['x2']:
        lam=(C['y2']-C['y1'])*pow(C['x2']-C['x1'],-1,P)%P
        X3=(lam*lam-C['x1']-C['x2'])%P; Y3=(lam*(C['x1']-X3)-C['y1'])%P
        print('   addition law: computed (x3,y3)=(%d,%d)'%(X3,Y3))
        print('   matches state x3,y3 ? %s / %s'%(X3==C['x3'],Y3==C['y3']))
        print('   offset x3_state - X3 = %d'%((C['x3']-X3)%P))
    D0=(v[7068]-v[2099])%P; K2=v[28730]%P
    print('   D0 = (x7068-x2099) mod p = %d'%D0)
    print('   K2 = x28730 mod p        = %d'%K2)
    print('   D0 is a valid secp256k1 x-coord (x^3+7 is QR): %s'%qr(pow(D0,3,P)+7))
    print('   K2 is a valid secp256k1 x-coord              : %s'%qr(pow(K2,3,P)+7))
    print('   b(D0,K2) = %d'%((K2*K2-pow(D0,3,P))%P))
    print('   b(D0,K2) == b1 ? %s ; == b2 ? %s ; == 7 ? %s'%(((K2*K2-pow(D0,3,P))%P)==b1,
          ((K2*K2-pow(D0,3,P))%P)==b2, ((K2*K2-pow(D0,3,P))%P)==7))
    for name,val in [('D0',D0),('K2',K2)]:
        for cn,cv in [('x1',C['x1']),('y1',C['y1']),('x2',C['x2']),('y2',C['y2']),('x3',C['x3']),('y3',C['y3'])]:
            if val==cv: print('   !! %s == %s'%(name,cn))
    print()
