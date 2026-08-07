import sys, os, pickle, collections, time, random
exec(open('/home/user/integer_solver/solve_lab/agentV_work/mirror/L/calib.py').read().split("if __name__=='__main__':")[0])
C=pickle.load(open('calib.pkl','rb')); ORIENT=C['ORIENT']
if __name__=='__main__':
    rnd=random.Random(11)
    for S in ([live[0]],[live[0],live[1]],rnd.sample(live,3),rnd.sample(live,17),rnd.sample(live,89),list(live)):
        v,isl,valn,order=assignment(set(S),ORIENT)
        vv,r=run(v)
        nz=[i for i,x in enumerate(r) if x%p]
        print('|S|=%-4d nonzero atoms mod p: %-4d  root=%s'%(len(S),len(nz),valn[ROOTID]))
        if 0<len(nz)<=25:
            for i in nz: print('      ',E.res[i][:150])
