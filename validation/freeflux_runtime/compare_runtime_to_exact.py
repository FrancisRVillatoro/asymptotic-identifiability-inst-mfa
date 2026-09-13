#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import csv
import json
import sys
import time
import warnings

import numpy as np
from scipy.linalg import expm

warnings.filterwarnings('ignore', category=FutureWarning)

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT))
from source_faithful_freeflux_runtime import build_toy_model, get_natural_MDV

TARGETS = ['Glu_123', 'Glu_12345', 'Cit_2345']
TARGET_TIMES = np.array([0.0,0.1,0.2,0.5,1.0,2.0])
P_NAT = 0.0107
POOLS = ['OAA','Cit','AKG','Suc','Fum','Glu']
NC = {'OAA':4,'Cit':6,'AKG':5,'Suc':4,'Fum':4,'Glu':5}
NISO = {m:2**NC[m] for m in POOLS}
offset={}; nstate=0
for m in POOLS:
    offset[m]=nstate; nstate += NISO[m]
CONST=nstate; N=nstate+1


def bits(i,n): return tuple((i>>(n-1-j))&1 for j in range(n))
def ibits(b):
    x=0
    for q in b: x=(x<<1)|int(q)
    return x

def trans_matrix(dst_n,src_n,mapper):
    T=np.zeros((2**dst_n,2**src_n))
    for j in range(2**src_n):
        for out,w in mapper(bits(j,src_n)):
            T[ibits(out),j]+=w
    return T

def bernoulli_iso(n,p=P_NAT):
    arr=np.zeros(2**n)
    for i in range(2**n):
        b=bits(i,n)
        arr[i]=np.prod([p if z else 1-p for z in b])
    return arr

# Exact positional representation of the FreeFlux AcCoA strategy.
_accoa=[]
for e in (0,1):
    pe=P_NAT if e else 1-P_NAT
    for f in (0,1):
        pf=P_NAT if f else 1-P_NAT
        _accoa.append(((e,f),0.5*pe*pf))
    _accoa.append(((e,1),0.25*pe))
_accoa.append(((1,1),0.25))
_acc={}
for b,w in _accoa: _acc[b]=_acc.get(b,0.0)+w
ACCOA_MIX=list(_acc.items())
assert abs(sum(w for _,w in ACCOA_MIX)-1)<1e-14

def map_v1(o):
    a,b,c,d=o
    return [((d,c,b,f,e,a),w) for (e,f),w in ACCOA_MIX]
T_v1=trans_matrix(6,4,map_v1)
T_v2=trans_matrix(5,6,lambda x:[(x[:5],1.)])
T_v3=np.eye(32)
def map_v4(x):
    s=x[1:5]
    return [(s,.5),(s[::-1],.5)]
T_v4=trans_matrix(4,5,map_v4)
def map_sym4(x): return [(x,.5),(x[::-1],.5)]
T_sym=trans_matrix(4,4,map_sym4)
EXT_OAA=bernoulli_iso(4)

def rs(m): return slice(offset[m],offset[m]+NISO[m])

def fluxes_from_theta(theta):
    u,w,r=np.exp(theta[:3])
    return {'v1':u+w,'v2':u+w,'v3':u,'v4':w,'v5':w,'v6_f':w+r,'v6_b':r,'v7':u}

def concentrations_from_theta(theta):
    return {m:np.exp(theta[3+i]) for i,m in enumerate(POOLS)}

def build_K(theta):
    F=fluxes_from_theta(theta); C=concentrations_from_theta(theta)
    K=np.zeros((N,N))
    def tr(dst,src,v,T): K[rs(dst),rs(src)] += (v/C[dst])*T
    def out(pool,v): K[rs(pool),rs(pool)] -= (v/C[pool])*np.eye(NISO[pool])
    def ext(dst,v,d): K[rs(dst),CONST] += (v/C[dst])*d
    tr('OAA','Fum',F['v6_f'],T_sym); ext('OAA',F['v7'],EXT_OAA); out('OAA',F['v1']+F['v6_b'])
    tr('Cit','OAA',F['v1'],T_v1); out('Cit',F['v2'])
    tr('AKG','Cit',F['v2'],T_v2); out('AKG',F['v3']+F['v4'])
    tr('Suc','AKG',F['v4'],T_v4); out('Suc',F['v5'])
    tr('Fum','Suc',F['v5'],T_sym); tr('Fum','OAA',F['v6_b'],T_sym); out('Fum',F['v6_f'])
    tr('Glu','AKG',F['v3'],T_v3); out('Glu',F['v3'])
    return K

def initial_state():
    q=np.zeros(N)
    for m in POOLS: q[rs(m)]=bernoulli_iso(NC[m])
    q[CONST]=1.
    return q
Q0=initial_state()

def mdv_matrix(pool,pos):
    pos=[p-1 for p in pos]
    M=np.zeros((len(pos)+1,N))
    for i in range(2**NC[pool]):
        b=bits(i,NC[pool]); M[sum(b[j] for j in pos),offset[pool]+i]+=1
    return M
OBS={
    'Glu_123':mdv_matrix('Glu',[1,2,3]),
    'Glu_12345':mdv_matrix('Glu',[1,2,3,4,5]),
    'Cit_2345':mdv_matrix('Cit',[2,3,4,5]),
}

theta0=np.array([np.log(5.),np.log(5.),np.log(7.5)]+[np.log(x) for x in [.1,5.,.3,1.,.2,.5]])
PARAM_NAMES=['log v3','log v4','log v6b']+[f'log C_{m}' for m in POOLS]

def exact_output(theta=theta0):
    K=build_K(theta)
    result={emu:{} for emu in TARGETS}
    for t in TARGET_TIMES:
        q=expm(K*t)@Q0
        for emu in TARGETS: result[emu][float(t)]=OBS[emu]@q
    return result

def make_grid(max_step):
    # Include every target exactly, subdividing each official interval.
    pts=[0.0]
    for a,b in zip(TARGET_TIMES[:-1],TARGET_TIMES[1:]):
        n=int(np.ceil((b-a)/max_step))
        pts.extend(np.linspace(a,b,n+1)[1:].tolist())
    return sorted(set(round(x,14) for x in pts))

def runtime_model():
    model,sim=build_toy_model(ROOT/'reactions.tsv')
    sim.prepare()
    return model,sim

def set_runtime_theta(model,theta):
    F=fluxes_from_theta(theta); C=concentrations_from_theta(theta)
    model.total_fluxes.loc[list(F)] = list(F.values())
    model.concentrations.loc[list(C)] = list(C.values())

def runtime_output(model,sim,theta,time_grid):
    set_runtime_theta(model,theta)
    model.timepoints=list(time_grid)
    raw=sim.simulate_raw()
    return {emu:{float(t):np.asarray(raw[emu][float(t)]) for t in TARGET_TIMES} for emu in TARGETS}

def stack(result,drop_last=False):
    vals=[]
    for emu in TARGETS:
        for t in TARGET_TIMES:
            x=np.asarray(result[emu][float(t)])
            vals.extend(x[:-1] if drop_last else x)
    return np.asarray(vals)

def compare(a,b):
    da=stack(a)-stack(b)
    return float(np.max(np.abs(da))),float(np.sqrt(np.mean(da**2))),float(np.linalg.norm(da)/np.linalg.norm(stack(b)))

def finite_difference(fun,theta,h=1e-5,drop_last=True):
    y0=stack(fun(theta),drop_last=drop_last)
    J=np.empty((len(y0),len(theta)))
    for j in range(len(theta)):
        tp=theta.copy(); tm=theta.copy(); tp[j]+=h; tm[j]-=h
        J[:,j]=(stack(fun(tp),drop_last=drop_last)-stack(fun(tm),drop_last=drop_last))/(2*h)
    return J

def main():
    exact=exact_output(theta0)
    model,sim=runtime_model()

    grids={
        'official':[0.0,0.1,0.2,0.5,1.0,2.0],
        'h0.1':make_grid(.1),
        'h0.05':make_grid(.05),
        'h0.025':make_grid(.025),
        'h0.0125':make_grid(.0125),
        'h0.00625':make_grid(.00625),
    }
    outputs={}
    conv_rows=[]
    for name,grid in grids.items():
        t0=time.time(); out=runtime_output(model,sim,theta0,grid); elapsed=time.time()-t0
        outputs[name]=out
        ma,rms,rel=compare(out,exact)
        conv_rows.append((name,len(grid),max(np.diff(grid)),ma,rms,rel,elapsed))
        print(name,len(grid),ma,rms,rel,elapsed)

    with open(OUT/'runtime_convergence.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['grid','n_timepoints','max_step','max_abs_vs_exact','rms_vs_exact','relative_Frobenius','seconds']); w.writerows(conv_rows)

    # Full target-by-target trajectory table for official grid and finest grid.
    with open(OUT/'trajectory_comparison.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['emu','time','mass','exact','runtime_official','runtime_refined','official_minus_exact','refined_minus_exact'])
        for emu in TARGETS:
            for t in TARGET_TIMES:
                for k,(e,o,r) in enumerate(zip(exact[emu][float(t)],outputs['official'][emu][float(t)],outputs['h0.00625'][emu][float(t)])):
                    w.writerow([emu,t,k,e,o,r,o-e,r-e])

    # Check official rounded fixture for Glu_12345.
    fixture={
        0.0:[.948,.051,.001,0,0,0], 0.1:[.928,.060,.012,0,0,0],
        0.2:[.873,.085,.041,.001,0,0], 0.5:[.697,.162,.131,.008,.001,0],
        1.0:[.518,.234,.213,.028,.006,0], 2.0:[.381,.274,.262,.063,.018,.001],
    }
    fixture_rows=[]; round_match=True
    for t in TARGET_TIMES:
        pred=outputs['official']['Glu_12345'][float(t)]
        expected=np.asarray(fixture[float(t)])
        match=np.array_equal(np.round(pred,3),expected)
        round_match &= match
        fixture_rows.append((t,*pred,*expected,float(np.max(np.abs(pred-expected))),match))
    with open(OUT/'official_fixture_check.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['time']+[f'runtime_m{i}' for i in range(6)]+[f'fixture_m{i}' for i in range(6)]+['max_abs','rounds_exactly']); w.writerows(fixture_rows)

    # Sensitivities: refined source-faithful runtime vs exact isotopomer solution.
    sens_step=.0125
    sens_grid=make_grid(sens_step)
    runtime_fun=lambda th: runtime_output(model,sim,th,sens_grid)
    exact_fun=lambda th: exact_output(th)
    print('finite-difference exact...')
    J_exact=finite_difference(exact_fun,theta0,h=1e-5,drop_last=True)
    print('finite-difference runtime...')
    J_run=finite_difference(runtime_fun,theta0,h=1e-5,drop_last=True)
    dJ=J_run-J_exact
    srows=[]
    for j,name in enumerate(PARAM_NAMES):
        absmax=float(np.max(np.abs(dJ[:,j])))
        rms=float(np.sqrt(np.mean(dJ[:,j]**2)))
        rel=float(np.linalg.norm(dJ[:,j])/max(np.linalg.norm(J_exact[:,j]),1e-300))
        srows.append((name,absmax,rms,rel,float(np.linalg.norm(J_exact[:,j])),float(np.linalg.norm(J_run[:,j]))))
    with open(OUT/'sensitivity_comparison.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['parameter','max_abs','rms','relative_L2','exact_L2','runtime_L2']); w.writerows(srows)
    np.savez_compressed(OUT/'sensitivity_matrices.npz',J_exact=J_exact,J_runtime=J_run,parameter_names=np.asarray(PARAM_NAMES))

    # Richardson estimate for the trajectory order.
    hs=np.array([row[2] for row in conv_rows[1:]],float)
    errs=np.array([row[3] for row in conv_rows[1:]],float)
    slope=float(np.polyfit(np.log(hs[-4:]),np.log(errs[-4:]),1)[0])

    summary={
        'freeflux_version': '0.3.8',
        'git_commit': 'ec05c47bbc2e4ac58bb85d39408ff4ef4016a15c',
        'official_fixture_round_match': bool(round_match),
        'official_grid_max_abs_vs_exact': conv_rows[0][3],
        'finest_grid_max_abs_vs_exact': conv_rows[-1][3],
        'finest_grid_rms_vs_exact': conv_rows[-1][4],
        'observed_convergence_order_max_error': slope,
        'sensitivity_internal_max_step': sens_step,
        'sensitivity_global_max_abs': float(np.max(np.abs(dJ))),
        'sensitivity_global_rms': float(np.sqrt(np.mean(dJ**2))),
        'sensitivity_global_relative_Frobenius': float(np.linalg.norm(dJ)/np.linalg.norm(J_exact)),
        'sensitivity_worst_parameter_relative_L2': max((r[3],r[0]) for r in srows),
        'target_emus': TARGETS,
        'target_times': TARGET_TIMES.tolist(),
        'compatibility_changes': ['np.float -> float','scipy.linalg.pinv2 -> pinv','fitting/optimization/display layers omitted'],
    }
    (OUT/'block1_summary.json').write_text(json.dumps(summary,indent=2))
    txt=['FreeFlux 0.3.8 source-faithful runtime validation','='*64]
    for k,v in summary.items(): txt.append(f'{k}: {v}')
    txt += ['', 'Per-parameter sensitivity comparison:']
    for r in srows: txt.append(f'  {r[0]:12s} max={r[1]:.3e} rms={r[2]:.3e} relL2={r[3]:.3e}')
    (OUT/'block1_summary.txt').write_text('\n'.join(txt))
    print('\n'.join(txt))

if __name__=='__main__':
    main()
