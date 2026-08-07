"""Export the curve/chain identification as JSON, with an independent re-verification."""
import os, sys, json, pickle
sys.path.insert(0,'/home/user/integer_solver/solve_lab/agentG_work')
import gpt
from gsym2 import L, ad, P
ch=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/chain.pkl','rb')); pts=pickle.load(open('/home/user/integer_solver/solve_lab/agentG_work/bitpoints.pkl','rb'))
order=ch['order']; P0=ch['P0']
Q=[None]*256; R=P0
for i in range(256): Q[i]=R; R=gpt.add(R,R)
bad=sum(1 for i,b in enumerate(order) if pts[b][0]!=Q[i])
import g46_table as T
base=T.frame([]); B1,B2,B3=base['pts']
out=dict(
 p=P, n=gpt.n,
 curve_in_instance="y^2 = x^3 + a2 x^2 + a4 x + a6",
 a2=gpt.K,
 a4=(gpt.K*gpt.K*pow(3,-1,P))%P,
 a6=77755683306591771556999954628254672912734268662742093169295805431582354953490,
 j_invariant=0,
 short_form_B=64019533680030876408443198762210829058751700634554282185987325820393598524794,
 iso_to_secp256k1="x_sec = (x + a2/3) / u^2 , y_sec = y / u^3",
 u=gpt._u,
 P0=list(P0), P1=list(B1), P2=list(B2), P3_target=list(B3),
 chain_bit_order=order,
 chain_verified_positions=256-bad,
 current_message_bits={"x2081":order.index(2081),"x24601":order.index(24601)},
 statement="message bit at chain index i pins the point [2^i]P0; the circuit adds the "
           "selected points and checks the sum equals P3.  A full solve requires "
           "k with [k]P0 = P3, i.e. the discrete logarithm of P3 to base P0 on secp256k1.")
json.dump(out,open('/home/user/integer_solver/solve_lab/agentG_work/secp_identification.json','w'),indent=1)
print('chain positions verified as [2^i]P0: %d/256'%(256-bad))
print('P1 == [2^%d]P0 : %s'%(order.index(24601), B1==Q[order.index(24601)]))
print('P2 == [2^%d]P0 : %s'%(order.index(2081),  B2==Q[order.index(2081)]))
print('P3 on curve, order n : %s , %s'%((B3[1]*B3[1]-pow(B3[0],3,P)-7)%P==0, gpt.mul(gpt.n,B3) is None))
print('written secp_identification.json')
