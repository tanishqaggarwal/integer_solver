import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
p=2**256-2**32-977
K1=33472904810391811973223207617762334363023286939839396241234196646906030803538671321618319   # x_6418 pin
K2=42775533402728869434716629464193396056515231264222641773817154079369026410240838606908039   # x_12553 pin
K3=119562606790549640390870952418684367882170154220603339634805704742270834564330392414192110  # x_31861 pin
K4=113141528427610260107049117992526537105383080782811760722361109500341947028737388716982706  # x_14865 pin
x24453=97553848499418123410591666447050222001188385549510401465815187079080512838891
# x_17925 = K1-K3, x_27019 = K2-K4  (constants)
d13 = K1-K3
d24 = K2-K4
A = d13*d13          # x_7010
B = d24*d24          # x_37530
C = K3+K1+x24453     # so x_7181 = x_9118 + C
# x_27177 = A*(x_9118+C) - B  ;  =0 => x_9118 = B/A - C
print('A=(K1-K3)^2, B=(K2-K4)^2')
print('B % A ==0 ?', B % A ==0)
print('gcd-based: does (K1-K3) | (K2-K4)?', (d24 % d13)==0, ' d24%d13=',d24%d13)
# General: x_27177=0 requires A | B. Check.
if B % A == 0:
    x9118 = B//A - C
    print('x_27177=0 solvable EXACTLY. x_9118 =', x9118)
else:
    print('x_27177=0 NOT exactly solvable (A does not divide B).')
    # fall back: we need x_31731=0 exactly and x_9106,x_2239 == 0 mod p.
# Let's also frame everything mod p and see if a joint solution exists.
# x_27177 = A*x_9118 + (A*C - B)   (linear in x_9118)
# x_4306  = d13*x_8731 + d24*x_9118 + (K4*d13 - d24*K3)
a27 = A % p ; b27 = (A*C - B) % p           # x_27177 ≡ a27*x_9118 + b27
# coefficients of x_4306: d13*x_8731 + d24*x_9118 + e ; e = K4*d13 - d24*K3
e4306 = (K4*d13 - d24*K3)
# loads mod p:
# x_31731 = 15964591*x_27177 + 13881285*x_4306  (want ==0 EXACT, but check mod p first)
# x_9106  = 7204959*x_27177  + 6822253*x_4306   (want ==0 mod p)
# x_2239  = 3494591*x_27177  + 14240157*x_4306  (want ==0 mod p)
print('--- checking joint mod-p feasibility ---')
# unknowns u=x_27177 mod p, w=x_4306 mod p. want:
#  9106:  7204959*u + 6822253*w ==0
#  2239:  3494591*u + 14240157*w ==0
# two homogeneous eqs => only u=w=0 unless determinant ==0 mod p
det = (7204959*14240157 - 6822253*3494591) % p
print('det of [9106;2239] system mod p =', det, ' -> only u=w=0' if det!=0 else ' dependent')
