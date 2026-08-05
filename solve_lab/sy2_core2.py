import os,sys; os.chdir('/home/user/integer_solver/solve_lab'); sys.path.insert(0,'.')
import heal_harness as H, json
p=H.p
fc=H.loadd('fc_partial.json')
for v in H.freeinp: H.val[v]=fc.get(v,0)
H.val[14865]=H.val[12553]; H.val[31861]=H.val[6418]; H.val[33168]=0
H.forward()
# sinks
print('x_37720=',H.val[37720],' /9994531 rem=',H.val[37720]%9994531)
print('x_9629=',H.val[9629],'x_23754=',H.val[23754],'x_35619=',H.val[35619])
H.val[950]=0; H.val[6947]=0
# x_8976 for atom 17901: x_37720 = 9994531*x_8976 -> x_8976 = x_37720/9994531 if divisible
if H.val[37720]%9994531==0:
    H.val[8976]=H.val[37720]//9994531
H.forward()
F=H.fails()
print('after zeroing sinks: fails=',len(F))
print(sorted(F))
# which atoms nonzero in remaining fails?
