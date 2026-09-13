#!/usr/bin/env python3
"""Run the unmodified FreeFlux 0.3.8 package in a supported Python 3.10 environment."""
from pathlib import Path
import csv
import numpy as np
from freeflux import Model, __version__

ROOT=Path(__file__).resolve().parent
assert __version__ == '0.3.8', __version__
model=Model('toy')
model.read_from_file(ROOT/'reactions.tsv')
isim=model.simulator('inst')
isim.set_target_EMUs({'Glu':[[1,2,3],'12345'],'Cit':'2345'})
isim.set_fluxes_from_file(ROOT/'fluxes.tsv')
isim.set_concentrations_from_file(ROOT/'concentrations.tsv')
isim.set_timepoints([0,0.1,0.2,0.5,1,2])
isim.set_labeling_strategy('AcCoA',labeling_pattern=['01','11'],percentage=[0.25,0.25],purity=[1,1],label_atom='C')
isim.prepare(n_jobs=1)
raw=isim.calculator._calculate_inst_MDVs()
rows=[]
for emu in ['Glu_123','Glu_12345','Cit_2345']:
    init=np.asarray(model.initial_sim_MDVs[emu][0])
    for k,x in enumerate(init): rows.append([emu,0,k,float(x)])
    for t in [0.1,0.2,0.5,1,2]:
        for k,x in enumerate(raw[emu][t]): rows.append([emu,t,k,float(x)])
with (ROOT/'official_freeflux_output.csv').open('w',newline='') as f:
    w=csv.writer(f); w.writerow(['emu','time','mass','value']); w.writerows(rows)
print('FreeFlux',__version__,'wrote official_freeflux_output.csv')
