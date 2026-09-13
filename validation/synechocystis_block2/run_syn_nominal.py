import sys, warnings, csv, json
from pathlib import Path
import openpyxl
warnings.filterwarnings('ignore', category=FutureWarning)
sys.path.insert(0,str(Path(__file__).resolve().parent))
import source_faithful_freeflux_runtime as ff

REA=str(Path(__file__).resolve().parent/'reactions.tsv')
CONC=str(Path(__file__).resolve().parent/'concentrations.xlsx')
fluxes = {
'pgi_f':0.13083510259427519,'pgi_b':3.7166874121552858,'g6pdh':3.2936325681668626,'pfk_f':4.5439874064737529,'pfk_b':12.788179945170489,'fba_f':0.1193413361379333,'fba_b':8.3635338748346317,'tpi_f':10.293426129319759,'tpi_b':21.670057167765151,'gapdh_f':3.2312406593978205,'gapdh_b':29.941353935507749,'gpm_f':30.39455603368928,'gpm_b':27.92450636549562,'eno_f':3.0886696473401001,'eno_b':0.6186199791464464,'pk_f':2.360913471579039,'pk_b':1.18835443501743,'rpe_f':2.0461647247035599,'rpe_b':9.836943453588102,'rpi_f':2.6696114770845929,'rpi_b':6.3281090192631613,'prk':14.742908839229941,'rbc1':14.702908839229869,'tkt1_f':2.7483219657733651,'tkt1_b':10.539100694657909,'tkt2_f':16.15530461871916,'tkt2_b':12.198453419500041,'tkt3_f':4.9055215207135667,'tkt3_b':1.071593991048164,'tal1_f':214.28199766210909,'tal1_b':213.58050863219242,'tal2_f':30.019590124833989,'tal2_b':30.721079154750743,'sba':3.1324384997485879,'sbp':3.1324384997485879,'pdh':1.3788671554123599,'cs':0.40610973149561841,'can_f':0.43199260173951892,'can_b':0.025882870243816653,'icd_f':0.67846214009517569,'icd_b':0.37593742041512601,'sdh_f':2.8563647254346201,'sdh_b':2.7527797136190038,'fum_f':2.0114817058436749,'fum_b':1.740318118568613,'mdh_f':15.059062578412769,'mdh_b':15.204313979322041,'icl':0.1035850118156447,'ms_f':7.8504397190700193,'ms_b':7.7668547072543461,'me':0.5000000000000403,'ppc':1.051643292524679,'rbc2':0.039999999999911162,'pgp':0.040000000000020956,'gld':0.040000000000020956,'gt':0.029999999999993438,'glyk_f':4.3045047601674655,'glyk_b':4.2745047601673809,'co2in':10.000000000000069,'biom':0.2453566258556989}
wb=openpyxl.load_workbook(CONC,data_only=True)
concs={r[0]:float(r[1]) for r in list(wb.active.iter_rows(values_only=True))[1:] if r[0] and r[1] is not None}
model=ff.Model('syn')
model.read_from_file(REA)
sim=ff.InstSimulator(model)
sim.set_target_EMUs({'G3P':['23','123'],'DHAP':'123','PEP':'123','Fum':'1234','R5P':'12345','RuBP':'12345','Cit':['12345','123456'],'S7P':'1234567'})
sim.set_labeling_strategy('CO2.ex', labeling_pattern=['1'], percentage=[0.5], purity=[0.997], label_atom='C')
for k,v in fluxes.items(): sim.set_flux(k,v)
for k,v in concs.items(): sim.set_concentration(k,v)
sim.set_timepoints([0,10,30,60,120,240])
print('prepare...')
sim.prepare()
print('EAM sizes', {k:len(v) for k,v in model.EAMs.items()})
print('simulate...')
out=sim.simulate_raw()
print('done targets',list(out))
# output CSV
rows=[]
for emu,d in out.items():
 for t,mdv in sorted(d.items()):
  for i,x in enumerate(mdv): rows.append((emu,t,i,float(x)))
with open(Path(__file__).resolve().parent/'nominal_mdv.csv','w',newline='') as f:
 w=csv.writer(f); w.writerow(['emu','time','mass','value']); w.writerows(rows)
print('n obs',len(rows))
for emu in out:
 print(emu, 't240', out[emu][240.0])
