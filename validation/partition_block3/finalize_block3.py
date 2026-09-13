#!/usr/bin/env python3
from pathlib import Path
import ast, importlib.util, json, glob, warnings
import numpy as np, pandas as pd, scipy.linalg as la, openpyxl
warnings.filterwarnings('ignore')
OUT=Path(__file__).resolve().parent; ROOT=OUT.parent/'synechocystis_block2'

# Load per-partition results
summ={}; vhs={}
for f in glob.glob(str(OUT/'block3_full_summary_*.json')):
    d=json.load(open(f)); summ.update(d['partitions'])
for f in glob.glob(str(OUT/'weak_subspaces_epsmin_*.npz')):
    z=np.load(f)
    for k in z.files: vhs[k]=z[k]

# Load benchmark/model for corrected nominal spectral diagnostics
spec=importlib.util.spec_from_file_location('ff',ROOT/'source_faithful_freeflux_runtime.py')
ff=importlib.util.module_from_spec(spec); spec.loader.exec_module(ff)
tr=ast.parse((ROOT/'run_syn_nominal.py').read_text()); fluxes=None
for node in tr.body:
    if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='fluxes' for t in node.targets): fluxes=ast.literal_eval(node.value); break
wb=openpyxl.load_workbook(ROOT/'concentrations.xlsx',data_only=True)
concs0={r[0]:float(r[1]) for r in list(wb.active.iter_rows(values_only=True))[1:] if r[0] and r[1] is not None}
model=ff.Model('syn'); model.read_from_file(ROOT/'reactions.tsv')
Stot=model._full_total_stoichiometric_matrix(tuple(model.metabolites),tuple(model.reactions)).drop(index=model.end_substrates+model.end_products)
total_ids=list(model.totalfluxids); f0=np.array([fluxes[k] for k in total_ids]); N=la.null_space(Stot.values); fss=N@(N.T@f0)
sim=ff.InstSimulator(model); sim.set_target_EMUs({'G3P':['23','123'],'DHAP':'123','PEP':'123','Fum':'1234','R5P':'12345','RuBP':'12345','Cit':['12345','123456'],'S7P':'1234567'}); sim.set_labeling_strategy('CO2.ex',labeling_pattern=['1'],percentage=[.5],purity=[.997],label_atom='C')
for k,v in zip(total_ids,fss): sim.set_flux(k,float(v))
for k,v in concs0.items(): sim.set_concentration(k,float(v))
sim.set_timepoints([0,10,30,60,120,240]); sim.prepare()
turn=pd.read_csv(OUT/'turnover_spectrum.csv'); turn_ser=turn.set_index('metabolite').turnover

def corrected_spectral(fast):
    rates=[]; details=[]
    for size in sorted(model.matrix_As):
        lambA,fa,pros=model.matrix_As[size]; lambM,mids=model.matrix_Ms[size]
        A=np.asarray(lambA(*model.total_fluxes[fa]),float)
        # IMPORTANT: nominal concentrations, independent of whichever epsilon was last used in sweeps
        M=np.asarray(lambM(*np.array([concs0[x] for x in mids])),float)
        F=la.pinv(M)@A
        idx=[i for i,e in enumerate(pros) if e.metabolite_id in fast]
        if idx:
            rr=-np.real(la.eigvals(F[np.ix_(idx,idx)])); pos=rr[rr>1e-10]
            if len(pos):
                rates.extend(pos.tolist()); details.append({'emu_size':int(size),'n_fast_emu':len(idx),'min_decay':float(pos.min()),'max_decay':float(pos.max())})
    slow=[m for m in turn_ser.index if m not in fast]
    slowmax=max(turn_ser[m] for m in slow) if slow else np.nan
    fastmin=min(turn_ser[m] for m in fast) if fast else np.nan
    return {'pool_gap_ratio':float(fastmin/slowmax) if slow else np.nan,
            'emu_fast_min_decay':float(min(rates)) if rates else np.nan,
            'emu_to_slow_turnover_ratio':float(min(rates)/slowmax) if rates and slow else np.nan,
            'spectral_by_size':details}

# Principal angles to primary biochemical partition
ref='auto_primary_phys'; sref=summ[ref]; idxref=[i-1 for i in sref['m1_indices_desc']]; Bref=vhs[ref][idxref].T
angle_rows=[]
for p,s in summ.items():
    idx=[i-1 for i in s['m1_indices_desc']]
    if idx:
        B=vhs[p][idx].T; a=np.degrees(la.subspace_angles(Bref,B)); amin=float(a.min()); amax=float(a.max()); comp=len(a)
    else: amin=amax=np.nan; comp=0
    s['angle_to_primary_phys_min_deg']=amin; s['angle_to_primary_phys_max_deg']=amax
    sd=corrected_spectral(s['fast_pools']); s.update(sd)
    angle_rows.append({'partition':p,'m1_dim':len(idx),'compared_dim':comp,'min_angle_deg':amin,'max_angle_deg':amax})

# Combined summary
rows=[]
for p,s in summ.items():
    other=s.get('other_indices_desc',[])
    # empirical robustness class
    if p in ('auto_primary_raw','auto_primary_phys','primary_plus1_raw','primary_plus1_phys','auto_fast1_raw') and len(other)==0 and len(s['m2_indices_desc'])==0:
        cls='quantized-stable'
    elif p.startswith('auto_secondary'):
        cls='broad-cut/crossover'
    else: cls='review'
    rows.append({'partition':p,'fast_pools':','.join(s['fast_pools']),'n_fast':s['n_fast'],
                 'm0_count':len(s['m0_indices_desc']),'m1_count':len(s['m1_indices_desc']),
                 'm1_exponents':';'.join(f'{x:.6f}' for x in s['m1_exponents']),
                 'm2_count':len(s['m2_indices_desc']),'other_count':len(other),'structural_count':s['structural_count'],
                 'pool_gap_ratio':s['pool_gap_ratio'],'emu_fast_min_decay':s['emu_fast_min_decay'],
                 'emu_to_slow_turnover_ratio':s['emu_to_slow_turnover_ratio'],
                 'output_convergence_order':s['output_convergence_order'],'output_last_relative_delta':s['output_last_relative_delta'],
                 'angle_to_primary_phys_max_deg':s['angle_to_primary_phys_max_deg'],'classification':cls})
pd.DataFrame(rows).sort_values(['classification','n_fast']).to_csv(OUT/'partition_stability_summary.csv',index=False)
pd.DataFrame(angle_rows).to_csv(OUT/'weak_subspace_angles.csv',index=False)

# Automatic gap candidates CSV from JSON
cands=json.load(open(OUT/'automatic_gap_candidates.json'))
pd.DataFrame(cands).to_csv(OUT/'automatic_gap_candidates.csv',index=False)

# Validate same gap rule on toy TCA benchmark.
toy_turn={'Cit':2.0,'Suc':5.0,'Glu':10.0,'AKG':100/3,'Fum':62.5,'OAA':175.0}
toy=pd.DataFrame({'metabolite':list(toy_turn),'turnover':list(toy_turn.values())}).sort_values('turnover').reset_index(drop=True)
toy['log10_turnover']=np.log10(toy.turnover); toy['gap_to_next']=toy.log10_turnover.shift(-1)-toy.log10_turnover; toy['ratio_to_next']=10**toy.gap_to_next
qmax=max(3,int(np.ceil(len(toy)/4))); tc=[]
for i in range(len(toy)-1):
    g=float(toy.loc[i,'gap_to_next']); prev=float(toy.loc[i-1,'gap_to_next']) if i else -np.inf; nxt=float(toy.loc[i+1,'gap_to_next']) if i<len(toy)-2 else -np.inf; suffix=len(toy)-i-1
    if suffix<=qmax and g>=prev and g>=nxt:
        tc.append({'cut_after':toy.loc[i,'metabolite'],'gap_log10':g,'ratio':10**g,'fast_size':suffix,'fast_pools':','.join(toy.loc[i+1:,'metabolite'])})
tc=sorted(tc,key=lambda d:d['gap_log10'],reverse=True); pd.DataFrame(tc).to_csv(OUT/'toy_gap_candidates.csv',index=False)
toy.to_csv(OUT/'toy_turnover_spectrum.csv',index=False)

# Machine-readable algorithm defaults
alg={'gap_metric':'g_i = log10(k_(i+1)/k_i)','tail_limit':'q_max = max(3, ceil(n/4))','candidate_rule':'local maxima of g_i among cuts leaving <= q_max pools in fast suffix','primary_rule':'largest such g_i','semantic_filter':'remove model-declared auxiliary/dummy atom-transfer pools before biochemical interpretation','spectral_audit':'eigenvalues of F_ff = M_ff^{-1} A_ff for each EMU size at nominal point','order_audit':'fit singular-value slopes over epsilon=1/64,1/128,1/256,1/512; classify m=0 if |slope|<0.15, m=1 if 0.75<slope<1.25, m=2 if 1.75<slope<2.25','robustness':'repeat primary raw, semantic-filtered, +/- boundary, and secondary local-gap candidates; compare principal angles of m=1 subspaces'}
(OUT/'partition_algorithm.json').write_text(json.dumps(alg,indent=2))

# Summary JSON
full={'algorithm':alg,'synechocystis_gap_candidates':cands,'toy_gap_candidates':tc,'partitions':summ}
(OUT/'BLOCK3_RESULTS.json').write_text(json.dumps(full,indent=2))
print(pd.DataFrame(rows).sort_values('n_fast').to_string(index=False))
print('\nangles')
print(pd.DataFrame(angle_rows).to_string(index=False))
print('\ntoy candidates')
print(pd.DataFrame(tc).to_string(index=False))
