"""Energy-limitation summary of param_sweep outputs: python analysis/energy_limitation_summary.py PARENT_DIR (one sub-folder per run set).
While substrate remains (> 5 %): days with an energy shortfall (R_endo), days with capped production, minimum E:C, share of uptake C
respired for energy, uptake-weighted CUE; plus day-365 % degraded, day of 50 % and maximum biomass."""
import pandas as pd, numpy as np, glob, sys
rows=[]
for d in sorted(glob.glob(sys.argv[1]+'/*')):
    lab=d.split('/')[-1]
    try: t=pd.read_csv(d+'/timeseries.csv.gz')
    except Exception: continue
    base=pd.read_csv(d+'/baseline_values.csv',header=None,index_col=0)[1]
    for p in base.index: t.loc[t.param=='baseline','value']=base[p]
    for (v,c,s),g in t.groupby(['value','climate','seed']):
        g=g.sort_values('doy'); dg=100-g.substrate_remaining_pct.values; dd=g.doy.values
        sub=g[g.substrate_remaining_pct>5]
        f=lambda x: dd[np.argmax(dg>=x)] if (dg>=x).any() else np.nan
        rows.append(dict(set=lab,value=v,climate=c.replace('_temp','').replace('_warm',''),deg=dg[-1],d50=f(50),
            def_days=int((sub.R_endo_C>0).sum()), capped=int(((sub.d_r_consti<0.999)|(sub.d_r_induci<0.999)).sum()),
            EC_min=sub.energy_to_C.min(), E_frac=1-np.average(sub.CAE.fillna(0),weights=sub.C_uptake.clip(lower=1e-12)),
            CUE=np.average(sub.CUE.fillna(0),weights=sub.C_uptake.clip(lower=1e-12)), bio_max=g.biomass_C.max(), bio_end=g.biomass_C.iloc[-1]))
m=pd.DataFrame(rows).groupby(['set','value','climate']).mean(numeric_only=True).unstack('climate')
pd.set_option('display.width',300); print(m.round(1).to_string())
