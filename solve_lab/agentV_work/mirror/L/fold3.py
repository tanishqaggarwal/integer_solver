import sys, pickle, random, collections
src=open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/calib2.py').read()
src=src.split("# numeric perm repair")[0]
exec(src)
C2=pickle.load(open('calib2.pkl','rb')); perm.update(C2['perm']); ORIENT=C2['ORIENT']
if __name__=='__main__':
    rnd=random.Random(3)
    for S in ([M['live'][0]],[M['live'][0],M['live'][1]],rnd.sample(M['live'],5),
              rnd.sample(M['live'],73),rnd.sample(M['live'],200),list(M['live'])):
        v,isl,valn=assignment(set(S),ORIENT)
        vv,r=run(v)
        nz=[i for i,x in enumerate(r) if x%p]
        print('|S|=%-4d nonzero atoms mod p: %-3d  root=%s'%(len(S),len(nz),valn[ROOT]))
        for i in nz: print('      ',E.res[i][:160])
