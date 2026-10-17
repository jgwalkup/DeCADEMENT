"""Grid summary for param_sweep outputs named p<production>_km<Vmax_Km> (one per sub-folder of PARENT_DIR), varying maintenance.
Per climate: day-365 % degraded, day 95 % degraded, energy-limited days while substrate > 5 % (and the first such day), minimum E:C,
biomass at day 365. Usage: python analysis/energy_limitation_grid_summary.py PARENT_DIR"""
import pandas as pd, numpy as np, glob, sys
rows=[]
for d in sorted(glob.glob(sys.argv[1]+'/p*')):
    try: t=pd.read_csv(d+'/timeseries.csv.gz')
    except Exception: continue
    b=d.split('/')[-1]; p=float(b[1:].split('_km')[0]); km=float(b.split('_km')[1])
    base=pd.read_csv(d+'/baseline_values.csv',header=None,index_col=0)[1].iloc[0]
    t.loc[t.param=='baseline','value']=base
    for (v,c,s),g in t.groupby(['value','climate','seed']):
        g=g.sort_values('doy'); dg=100-g.substrate_remaining_pct.values; dd=g.doy.values; sub=g[g.substrate_remaining_pct>5]
        lim=(sub.R_endo_C>0)|(sub.d_r_consti<0.999)|(sub.d_r_induci<0.999)
        rows.append(dict(prod=p,km=km,maint=v,clim=c.split('_')[0][:3],deg=dg[-1],exh=dd[np.argmax(dg>=95)] if (dg>=95).any() else np.nan,
            lim=int(lim.sum()),first=sub.doy[lim].min() if lim.any() else np.nan,ECmin=sub.energy_to_C.min() if len(sub) else np.nan,bio_end=g.biomass_C.iloc[-1]))
m=pd.DataFrame(rows).groupby(['prod','km','maint','clim']).mean(numeric_only=True).drop(columns='seed',errors='ignore').unstack('clim')
pd.set_option('display.width',260); print(m[['deg','exh','lim','first','ECmin','bio_end']].round(1).to_string())
