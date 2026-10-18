"""Global gain search: differential evolution + Nelder-Mead polish on the robust design cost.

   python3 tune_global.py base final_base.json    # PID, S-PID, FS-PID, QPID, S-QPID
   python3 tune_global.py eq   final_eq.json      # ES-QPID and its three ablations
   python3 lock_gains.py                          # round to 5 significant figures -> locked.json
Requires numpy, scipy, numba. Several hours on two cores.
"""
import json, sys, numpy as np
from multiprocessing import Pool
from scipy.optimize import differential_evolution, minimize
from sqpid_sim_py import *
SCEN=[dict(),dict(Lerr=0.9),dict(Lerr=1.1),dict(kg=0.9),dict(kg=1.1),dict(gamp=1.0),dict(t_on=3.0,t_off=9.0)]
W_OS=0.5; W_U=1.0; W_SS=5.0
B={'Kp':(-2,1.7),'Ki':(-3,2.7),'Kd':(-3,1),'Ks':(-2,1.7),'Tf':(-1.7,0.5),'mmax':(-2.3,0),'alpha':(-2,0),'wq':(-1.5,0.3),'npow':(0,1.2),'tm':(-1.3,0.69897)}
ORDER=['Kp','Ki','Kd','Ks','Tf']
def build(name,vals):
    law,sp,filt=CTRL[name]
    p=[vals['Kp'],vals['Ki'],vals['Kd']]
    if law>0: p.append(vals['Ks'])
    if filt: p.append(vals['Tf'])
    return p
def evalc(name,L,vals,opts,scen=SCEN,**extra):
    p=build(name,vals); o2=dict(opts or {})
    for k in ('mmax','alpha','wq','npow','tm'):
        if k in vals: o2[k]=vals[k]
    return [run(name,L,p,opts=o2,**s,**extra) for s in scen]
def cost(lp,name,L,tuned,fixed,opts):
    lp=np.asarray(lp); pen=0.0
    for i,k in enumerate(tuned):
        lo,hi=B[k]
        if lp[i]<lo: pen+=10*(lo-lp[i]); 
        if lp[i]>hi: pen+=10*(lp[i]-hi)
    lp=np.clip(lp,[B[k][0] for k in tuned],[B[k][1] for k in tuned])
    vals=dict(fixed); vals.update({k:10**lp[i] for i,k in enumerate(tuned)})
    J=0
    for o in evalc(name,L,vals,opts):
        if not np.isfinite(o['IAE']): return 1e3
        J+=o['sITAE']+W_OS*max(o['OS']-2,0)+W_U*o['su']+W_SS*max(o['ess']-0.02,0)+20*o['div']
    return J/len(SCEN)+pen
def job(a):
    tag,name,L,tuned,fixed,opts=a
    bnds=[B[k] for k in tuned]
    r=differential_evolution(cost,bnds,args=(name,L,tuned,fixed,opts),seed=1,popsize=16,maxiter=200,tol=1e-7,polish=False,init='sobol')
    r2=minimize(cost,r.x,args=(name,L,tuned,fixed,opts),method='Nelder-Mead',options=dict(maxiter=800,xatol=1e-4,fatol=1e-7))
    x=r2.x if r2.fun<r.fun else r.x
    x=np.clip(x,[b[0] for b in bnds],[b[1] for b in bnds])
    vals=dict(fixed); vals.update({k:float(10**x[i]) for i,k in enumerate(tuned)})
    return tag,name,L,vals,float(min(r.fun,r2.fun)),opts
def runjobs(jobs,out):
    res={}
    with Pool(2) as pl:
        for tag,name,L,vals,J,opts in pl.imap_unordered(job,jobs):
            res[f'{tag}|{L}']=dict(name=name,vals=vals,J=J,opts=opts)
            print(tag,L,{k:round(v,4) for k,v in vals.items()},'J=%.4f'%J,flush=True)
    json.dump(res,open(out,'w'),indent=1)
if __name__=='__main__':
    which=sys.argv[1]; Ls=(0.3,0.8,1.5); jobs=[]
    if which=='base':
        for L in Ls:
            jobs+=[('PID','PID',L,['Kp','Ki','Kd'],{},None),
                   ('S-PID','S-PID',L,['Kp','Ki','Kd'],{},None),
                   ('QPID','QPID',L,['Kp','Ki','Kd','Ks'],{},None),
                   ('S-QPID','S-QPID',L,['Kp','Ki','Kd','Ks'],{},None),
                   ('FS-PID','FS-PID',L,['Kp','Ki','Kd','Tf'],{},None)]
    elif which=='eq':
        for L in Ls:
            jobs+=[('ES-QPID','FS-QPID',L,['Kp','Ki','Kd','Ks','Tf','tm'],{},dict(ew=1,dwt=1)),
                   ('ES-noDw','FS-QPID',L,['Kp','Ki','Kd','Ks','Tf','tm'],{},dict(ew=1)),
                   ('ES-nofilt','S-QPID',L,['Kp','Ki','Kd','Ks','tm'],{},dict(ew=1,dwt=1)),
                   ('ES-noweight','FS-QPID',L,['Kp','Ki','Kd','Ks','Tf'],{},None)]
    runjobs(jobs,sys.argv[2])
