#!/usr/bin/env python3
"""Q-1: invariant of the stage law.
The law is  out_x = l^2 - a_x - b_x - K,  out_y = l*(a_x-out_x) - a_y,  l=(b_y-a_y)/(b_x-a_x).
Substituting X = x + c with 3c = K removes the offset:  out_X = l^2 - A_X - B_X.
That normalised law is the classical chord construction, whose CONSERVED SET is a cubic
   Y^2 = X^3 + aX + b.
Test, purely computationally: do the 256 leaf pin pairs lie on one such cubic?
"""
import json, os, itertools, sys
HERE = os.path.dirname(os.path.abspath(__file__))
p = 115792089237316195423570985008687907853269984665640564039457584007908834671663
K = 97553848499418123410591666447050222001188385549510401465815187079080512838891
c = K * pow(3, p-2, p) % p
pins = json.load(open(os.path.join(HERE, '..', 'agentF_work', 'pins.json')))
print('pins:', len(pins), 'shift c = K/3 =', c)

pts = []
for g, vv in pins.items():
    if len(vv) != 2: continue
    (v1, k1), (v2, k2) = vv
    pts.append((g, v1, k1 % p, v2, k2 % p))

def fit(P, Q):
    """a,b with Y^2 = X^3+aX+b through P,Q; None if X_P==X_Q."""
    (x1,y1),(x2,y2) = P,Q
    if (x1-x2) % p == 0: return None
    a = ((y1*y1-x1**3) - (y2*y2-x2**3)) % p * pow((x1-x2) % p, p-2, p) % p
    b = (y1*y1 - x1**3 - a*x1) % p
    return (a,b)

def on(a,b,X,Y): return (Y*Y - X**3 - a*X - b) % p == 0

# orientation: for each pin, is (k1,k2) = (x,y) or (y,x)?  Try the two global conventions.
for orient in (0,1):
    XY = []
    for g,v1,k1,v2,k2 in pts:
        X,Y = (k1,k2) if orient==0 else (k2,k1)
        XY.append((g,(X + c) % p, Y % p))
    ab = fit((XY[0][1],XY[0][2]), (XY[1][1],XY[1][2]))
    if ab is None: print('orient',orient,'degenerate'); continue
    a,b = ab
    good = sum(1 for g,X,Y in XY if on(a,b,X,Y))
    print('orient %d: a=%d\n           b=%d\n           on-curve %d/%d' % (orient,a,b,good,len(XY)))
    if good == len(XY):
        disc = (4*pow(a,3,p) + 27*b*b) % p
        print('  *** ALL LEAVES LIE ON ONE CUBIC ***')
        print('  discriminant 4a^3+27b^2 mod p =', disc)
        json.dump({'a':str(a),'b':str(b),'c_shift':str(c),'orient':orient,'p':str(p)},
                  open(os.path.join(HERE,'curve.json'),'w'))
