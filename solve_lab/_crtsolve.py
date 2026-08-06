import json
p=2**256-2**32-977
G=json.load(open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/gadget_affine.json'))
b={int(k):int(v) for k,v in G['b'].items()}
a9118={int(k):int(v) for k,v in G['s9118'].items()}
a8731={int(k):int(v) for k,v in G['s8731'].items()}
x7068=int(G['x7068']); x4432=int(G['x4432'])
# affine forms: x_q = b[q] + a9118[q]*x_9118 + a8731[q]*x_8731
# parametrize x_9118 = x7068 - 7376877*p*s ; x_8731 = x4432 - p*t   (s=x_17325, t=x_9413)
M1=7376877*p; M2=p
def form_in_st(q):
    # x_q = b[q] + a9118[q]*(x7068 - M1*s) + a8731[q]*(x4432 - M2*t)
    const = b[q] + a9118[q]*x7068 + a8731[q]*x4432
    cs = -a9118[q]*M1
    ct = -a8731[q]*M2
    return const, cs, ct   # x_q = const + cs*s + ct*t
c31,s31,t31 = form_in_st(31731)
c22,s22,t22 = form_in_st(2239)
c91,s91,t91 = form_in_st(9106)
print("x_31731(s,t) = %d + %d*s + %d*t"%(c31,s31,t31))
print("  need = 0 exactly")
print("x_2239 mod p: const=%d, cs=%d, ct=%d (mod p)"%(c22%p, s22%p, t22%p))
print("x_9106 mod (13523997p): checking...")
from math import gcd
# Constraint 1: x_31731 = 0 exactly: c31 + s31*s + t31*t = 0
g=gcd(s31,t31)
print("gcd(s31,t31)=",g," divides c31?", c31%g==0)
# Constraint 2: x_2239 ≡ 0 mod p
print("x_2239 const mod p =", c22%p, " s-coef mod p=", s22%p," t-coef mod p=", t22%p)
# Constraint 3: x_9106 ≡ 0 mod 13523997*p
M3=13523997*p
print("x_9106 const mod M3=", c91%M3, " s-coef mod M3=", s91%M3," t-coef mod M3=", t91%M3)
# Try to solve: from x_31731=0 exact, express one var. If s31 or t31 divides nicely.
# General: s31*s + t31*t = -c31. Solve with extended euclid.
def ext(a,b):
    if b==0: return (a,1,0)
    g,x,y=ext(b,a%b); return (g,y,x-(a//b)*y)
g2,u,v=ext(s31,t31)
if (-c31)%g2==0:
    s0=u*(-c31//g2); t0=v*(-c31//g2)
    # general sol: s=s0 + (t31/g2)*k, t=t0 - (s31/g2)*k
    ds=t31//g2; dt=-s31//g2
    print(f"\nx_31731=0 solution family: s={s0}+{ds}*k, t={t0}+{dt}*k")
    # substitute into constraint2 (mod p) and constraint3 (mod M3) -> congruences in k
    # x_2239 ≡ 0 mod p: c22 + s22*(s0+ds*k) + t22*(t0+dt*k) ≡ 0 mod p
    A2=(s22*ds + t22*dt)%p
    B2=(c22 + s22*s0 + t22*t0)%p
    print(f"constraint2: {A2}*k + {B2} ≡ 0 mod p")
    # x_9106 ≡ 0 mod M3
    A3=(s91*ds + t91*dt)%M3
    B3=(c91 + s91*s0 + t91*t0)%M3
    print(f"constraint3: {A3}*k + {B3} ≡ 0 mod M3")
    # solve constraint2 for k mod p
    sols=[]
    if A2%p!=0:
        k2=(-B2*pow(A2,p-2,p))%p
        print(f"k ≡ {k2} mod p (from constraint2)")
        # solve constraint3: A3*k ≡ -B3 mod M3
        g3=gcd(A3,M3)
        if (-B3)%g3==0:
            M3r=M3//g3
            k3=((-B3//g3)*pow((A3//g3)%M3r,-1,M3r))%M3r
            print(f"k ≡ {k3} mod {M3r} (from constraint3)")
            # CRT combine k≡k2 mod p, k≡k3 mod M3r
            from math import gcd as _g
            gg=_g(p,M3r)
            if (k2-k3)%gg==0:
                lcm=p//gg*M3r
                # CRT
                k=(k2 + p*(((k3-k2)//gg*pow((p//gg)%(M3r//gg),-1,M3r//gg))%(M3r//gg)))%lcm
                print(f"CRT: k ≡ {k} mod {lcm}  -> SOLUTION EXISTS")
                json.dump({'s0':s0,'t0':t0,'ds':ds,'dt':dt,'k':int(k),'lcm':int(lcm)},
                          open('/tmp/claude-0/-home-user-integer-solver/d279c5f7-7ad9-5092-a944-cb7784834ec0/scratchpad/crtsol.json','w'))
            else:
                print("CRT INCOMPATIBLE: (k2-k3) not divisible by gcd(p,M3r)=",gg, " -> NO SOLUTION (this is the wall)")
        else:
            print("constraint3 unsolvable: gcd doesn't divide")
    else:
        print("A2≡0 mod p; constraint2 is", B2%p==0 and "auto-satisfied" or "unsatisfiable")
else:
    print("x_31731=0 has no integer solution (gcd)")
