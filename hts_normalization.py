import argparse, os, re
from datetime import datetime
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm

def normalize_well(x):
    s=str(x).strip().upper(); m=re.fullmatch(r"([A-P])(\d{1,2})",s)
    return f"{m.group(1)}{int(m.group(2)):02d}" if m else s

def resolve_file(prompt, default=None):
    while True:
        v=input(f"{prompt} [{default}]: ").strip().strip('"') if default else input(f"{prompt}: ").strip().strip('"')
        v=v or default; p=Path(v).expanduser()
        if not p.exists() and p.suffix=='':
            q=Path(str(p)+'.csv')
            if q.exists(): p=q
        if p.is_file(): return p
        print(f"\nERROR: File not found: {p}\nEnter the filename including .csv, or the full path.\n")

def choose_column(df,prompt,default=None,numeric_only=False):
    cols=list(df.columns)
    candidates=[c for c in cols if pd.to_numeric(df[c],errors='coerce').notna().any()] if numeric_only else cols
    if not candidates: candidates=cols
    print("\nAvailable columns:"); [print(f"  {i}. {c}") for i,c in enumerate(candidates,1)]
    while True:
        v=input(f"{prompt} [{default}]: " if default is not None else f"{prompt}: ").strip() or (str(default) if default is not None else '')
        try:
            n=int(v)
            if 1<=n<=len(candidates): return candidates[n-1]
        except ValueError: pass
        print("Please enter a valid column number.")

def choose_optional(df,prompt):
    print("\nAvailable columns:"); [print(f"  {i}. {c}") for i,c in enumerate(df.columns,1)]; print("  0. None / not used")
    while True:
        try:
            n=int(input(f"{prompt} [0]: ").strip() or '0')
            if n==0:return None
            if 1<=n<=len(df.columns):return df.columns[n-1]
        except ValueError: pass
        print("Please enter a valid number.")

def controls(df,well_col,layout,pos,neg):
    s=pd.Series('',index=df.index,dtype='object')
    if layout=='columns':
        c=pd.to_numeric(df['_WellColumn'],errors='coerce'); ns={int(x.strip()) for x in neg.split(',') if x.strip()}; ps={int(x.strip()) for x in pos.split(',') if x.strip()}
        s.loc[c.isin(ns)]='NEGATIVE'; s.loc[c.isin(ps)]='POSITIVE'
    elif layout=='rows':
        r=df[well_col].str[0].str.upper(); s.loc[r.isin({x.strip().upper() for x in neg.split(',') if x.strip()})]='NEGATIVE'; s.loc[r.isin({x.strip().upper() for x in pos.split(',') if x.strip()})]='POSITIVE'
    else:
        w=df[well_col].map(normalize_well); s.loc[w.isin({normalize_well(x) for x in neg.split(',') if x.strip()})]='NEGATIVE'; s.loc[w.isin({normalize_well(x) for x in pos.split(',') if x.strip()})]='POSITIVE'
    return s

def plate_calc(g,assay):
    g=g.copy(); neg=g.loc[g['Control State']=='NEGATIVE',assay]; pos=g.loc[g['Control State']=='POSITIVE',assay]; nm,pm=neg.mean(),pos.mean(); sd=g[assay].std(); g['Raw data z-score']=(g[assay]-g[assay].mean())/sd if pd.notna(sd) and sd else None; den=nm-pm; g['Percent inhibition']=100-((g[assay]-pm)/den*100) if pd.notna(den) and den else None; isd=g['Percent inhibition'].std(); g['Percent inhibition z-score']=(g['Percent inhibition']-g['Percent inhibition'].mean())/isd if pd.notna(isd) and isd else None; return g

def qc(g,assay):
    n=g.loc[g['Control State']=='NEGATIVE',assay]; p=g.loc[g['Control State']=='POSITIVE',assay]
    if len(n)==0 or len(p)==0:return None,'UNKNOWN'
    den=abs(n.mean()-p.mean())
    if pd.isna(den) or den==0:return None,'UNKNOWN'
    z=1-3*(n.std()+p.std())/den; return z,('GOOD' if z>=.5 else ('BORDERLINE' if z>=0 else 'BAD'))

def analyze(a):
    d=pd.read_csv(a.data); b=pd.read_csv(a.database); d.columns=d.columns.astype(str).str.strip(); b.columns=b.columns.astype(str).str.strip()
    d[a.data_plate]=d[a.data_plate].astype(str).str.strip().str.upper(); d[a.data_well]=d[a.data_well].map(normalize_well); b[a.db_plate]=b[a.db_plate].astype(str).str.strip().str.upper(); b[a.db_well]=b[a.db_well].map(normalize_well); d[a.assay]=pd.to_numeric(d[a.assay],errors='coerce'); d['_WellColumn']=d[a.data_well].str[1:].astype(int)
    d['Control State']=controls(d,a.data_well,a.control_layout,a.positive_spec,a.negative_spec)
    m=d.merge(b,left_on=[a.data_plate,a.data_well],right_on=[a.db_plate,a.db_well],how='left',validate='many_to_one',suffixes=('','_database'))
    missing=m[a.db_compound].isna().sum(); print(f"\nWarning: {missing} wells did not match a compound record.") if missing else None
    f=pd.concat([plate_calc(g,a.assay) for _,g in m.groupby(a.data_plate,sort=False)],ignore_index=True); nc=f[f['Control State']=='']; mean=nc['Percent inhibition'].mean(); sd=nc['Percent inhibition'].std(); f['Run percent inhibition z-score']=(f['Percent inhibition']-mean)/sd if pd.notna(sd) and sd else None; f['Hit']=(f['Control State']=='')&(f['Percent inhibition']>=a.hit_threshold)
    q=pd.DataFrame([{'Plate':p,'Z-prime':qc(g,a.assay)[0],'QC':qc(g,a.assay)[1]} for p,g in f.groupby(a.data_plate,sort=False)])
    cols=[a.data_plate,a.data_well,a.db_compound,a.db_structure,a.assay,'Raw data z-score','Percent inhibition','Percent inhibition z-score','Run percent inhibition z-score','Hit','Control State']
    for c in (a.db_name,a.db_batch):
        if c and c in f.columns and c not in cols: cols.insert(2,c)
    out=f[[c for c in cols if c in f.columns]]; od=Path(a.output); od.mkdir(parents=True,exist_ok=True); date=datetime.fromtimestamp(os.path.getmtime(a.data)).strftime('%Y-%m-%d'); safe=re.sub(r'[^\w.-]+','_',a.assay).strip('_') or 'assay'; csv=od/f'{safe}_normalized_{date}.csv'; par=od/f'{safe}_normalized_{date}.parquet'; xlsx=od/f'{safe}_top_hits_{date}.xlsx'; qcsv=od/f'{safe}_plate_QC_{date}.csv'; pdf=od/f'{safe}_plate_report_{date}.pdf'; out.to_csv(csv,index=False); q.to_csv(qcsv,index=False)
    try: out.to_parquet(par,index=False)
    except Exception as e: print(f'Warning: Parquet not created ({e}). Install pyarrow for Parquet output.')
    with pd.ExcelWriter(xlsx,engine='openpyxl') as w: out[out['Hit']==True].sort_values('Percent inhibition',ascending=False).to_excel(w,sheet_name='Top Hits',index=False); q.to_excel(w,sheet_name='Plate QC',index=False)
    with PdfPages(pdf) as p:
        for plate,g in f.groupby(a.data_plate,sort=False):
            g=g.copy(); g['_Row']=g[a.data_well].str[0]; g['_Col']=g[a.data_well].str[1:].astype(int); rows=list('ABCDEFGHIJKLMNOP'); cs=list(range(1,25)); h=g.pivot(index='_Row',columns='_Col',values=a.assay).reindex(index=rows,columns=cs); i=g.pivot(index='_Row',columns='_Col',values='Percent inhibition').reindex(index=rows,columns=cs); z,_=qc(g,a.assay); fig,ax=plt.subplots(1,2,figsize=(16,8)); ax[0].imshow(h,aspect='auto',cmap='Blues_r'); ax[0].set_title(a.assay); ax[1].imshow(i,aspect='auto',cmap='coolwarm',norm=TwoSlopeNorm(vmin=-100,vcenter=0,vmax=100)); ax[1].set_title('Percent inhibition'); [x.set_xticks(range(24)) or x.set_xticklabels(cs) for x in ax]; [x.set_yticks(range(16)) or x.set_yticklabels(rows) for x in ax]; fig.suptitle(f'Plate {plate} — Z-prime = {z:.3f}' if z is not None else f'Plate {plate} — Z-prime = NA'); plt.tight_layout(rect=[0,0,1,.92]); p.savefig(fig); plt.close(fig)
    print(f"\nAnalysis complete.\nPlates analyzed: {len(q)}\nPotential hits: {int(f['Hit'].sum())}\nOutput folder: {od.resolve()}")

def interactive():
    print('\nDDRC High-Throughput Screening Normalization\n============================================')
    data=resolve_file('Screening dataset CSV','testdata.csv'); db=resolve_file('Compound database CSV','testdb.csv'); dp=pd.read_csv(data,nrows=20); cp=pd.read_csv(db,nrows=20); dp.columns=dp.columns.astype(str).str.strip(); cp.columns=cp.columns.astype(str).str.strip()
    print('\nSCREENING DATASET'); data_plate=choose_column(dp,'Which column contains the PLATE identifier?',1); data_well=choose_column(dp,'Which column contains the WELL identifier?',2); assay=choose_column(dp,'Which column contains the RAW DATA / READOUT measurement?',4,True); choose_optional(dp,'Which column contains existing control information? (optional)')
    print('\nCOMPOUND DATABASE'); db_plate=choose_column(cp,'Which column contains the database PLATE identifier?',4); db_well=choose_column(cp,'Which column contains the database WELL identifier?',5); db_compound=choose_column(cp,'Which column contains the COMPOUND identifier?',1); db_structure=choose_column(cp,'Which column contains the STRUCTURE / SMILES / CXSMILES?',2); db_name=choose_optional(cp,'Optional: which column contains a molecule/name description?'); db_batch=choose_optional(cp,'Optional: which column contains batch information?')
    print('\nCONTROL LAYOUT\n  1. Entire column(s)\n  2. Entire row(s)\n  3. Specific wells'); layout={'1':'columns','2':'rows','3':'wells'}.get(input('How are the positive and negative controls arranged? [1]: ').strip() or '1');
    if not layout: raise ValueError('Invalid control layout.')
    if layout=='columns': pos=input('Which column(s) contain the POSITIVE control? [e.g. 24]: ').strip(); neg=input('Which column(s) contain the NEGATIVE control? [e.g. 23]: ').strip()
    elif layout=='rows': pos=input('Which row(s) contain the POSITIVE control? [e.g. B]: ').strip(); neg=input('Which row(s) contain the NEGATIVE control? [e.g. A]: ').strip()
    else: pos=input('Which well(s) contain the POSITIVE control? [e.g. A01,A02]: ').strip(); neg=input('Which well(s) contain the NEGATIVE control? [e.g. B01,B02]: ').strip()
    hit=float(input('Hit threshold (% inhibition) [50]: ').strip() or '50'); output=input('Output directory [Output]: ').strip() or 'Output'; print(f'\nConfiguration: Plate={data_plate}; Well={data_well}; Raw data={assay}; DB Plate={db_plate}; DB Well={db_well}; Compound={db_compound}; Structure={db_structure}; Controls={layout}; Positive={pos}; Negative={neg}; Hit threshold={hit}%; Output={output}')
    if input('Continue? [Y/n]: ').strip().lower()=='n': raise SystemExit('Analysis cancelled.')
    return argparse.Namespace(data=data,database=db,data_plate=data_plate,data_well=data_well,assay=assay,db_plate=db_plate,db_well=db_well,db_compound=db_compound,db_structure=db_structure,db_name=db_name,db_batch=db_batch,control_layout=layout,positive_spec=pos,negative_spec=neg,hit_threshold=hit,output=output)

def main():
    p=argparse.ArgumentParser(); p.add_argument('--data'); p.add_argument('--database'); p.add_argument('--data-plate'); p.add_argument('--data-well'); p.add_argument('--assay'); p.add_argument('--db-plate'); p.add_argument('--db-well'); p.add_argument('--db-compound'); p.add_argument('--db-structure'); p.add_argument('--db-name'); p.add_argument('--db-batch'); p.add_argument('--control-layout',choices=['columns','rows','wells']); p.add_argument('--positive-spec'); p.add_argument('--negative-spec'); p.add_argument('--hit-threshold',type=float,default=50); p.add_argument('--output',default='Output'); a=p.parse_args(); analyze(a) if all([a.data,a.database,a.data_plate,a.data_well,a.assay,a.db_plate,a.db_well,a.db_compound,a.db_structure,a.control_layout,a.positive_spec,a.negative_spec]) else analyze(interactive())
if __name__=='__main__': main()
