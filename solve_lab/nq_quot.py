import heal_harness as H
p=H.p
d=H.loadd('sy_regime11_39018.json')
for v in H.freeinp: H.val[v]=d.get(v,0)
for k,vv in d.items(): H.val[k]=vv
H.forward()
V=H.val
def g(n): return V[n]
# gadget intermediate values
x27019=g(12553)-g(14865); x17925=g(6418)-g(31861)
x34310=g(31861)-g(9118); x3349=g(8731)+g(14865)
print(f"x_21279 (selector) = {g(21279)}")
print(f"x_27019 = x_12553-x_14865 = {x27019}  (|{len(str(abs(x27019)))}d|)")
print(f"x_17925 = x_6418-x_31861 = {x17925}  (|{len(str(abs(x17925)))}d|)")
print(f"x_34310 = x_31861-x_9118 = {x34310}")
print(f"x_3349  = x_8731+x_14865 = {x3349}")
print(f"x_7181 = {g(7181)}, x_27177={g(27177)}, x_4306={g(4306)}")
print(f"x_31731={g(31731)}, x_9106={g(9106)}")
# status of knobs
for n in [7181,8731,9118,6947,950,33168,2239,26874]:
    print(f"  x_{n}: {'FREE' if n in H.freeinp else 'gate'} = {g(n)}")
# quotient fix feasibility:
print("\n=== quotient fix divisibility checks ===")
# x_27177 = x_7181*x_17925^2 - x_27019^2 = 0 -> x_7181 = x_27019^2/x_17925^2
num=x27019*x27019; den=x17925*x17925
print(f"x_27019^2 % x_17925^2 == 0 ? {num%den==0}  -> x_7181 would be {num//den if den and num%den==0 else 'N/A'}")
# x_4306 = x_3349*x_17925 - x_27019*x_34310 = 0 -> x_34310 = x_3349*x_17925/x_27019 (via x_9118)
num2=x3349*x17925
print(f"x_3349*x_17925 % x_27019 == 0 ? {num2%x27019==0}  -> x_34310 target {num2//x27019 if num2%x27019==0 else 'N/A'}")
# check gcd relationships
import math
print(f"gcd(x_17925, x_27019) = {math.gcd(abs(x17925),abs(x27019))}")
