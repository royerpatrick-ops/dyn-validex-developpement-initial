import numpy as np, pandas as pd  # pyright: ignore[reportMissingImports]
d=pd.read_csv('balayage_translations.csv'); inv=pd.read_csv('invariants.csv')
d=d.merge(inv,on=['forme','param'])
d['C']=d.deficit_abs*d.R**(1/3)/d.I43
d['cm']=d.m_sommets/(d.Omega*d.R**(2/3))
g=d.groupby(['forme','param','R']).agg(C=('C','mean'),Cse=('C','sem'),cm=('cm','mean')).reset_index()
lab=lambda r: f"{'super p=' if r.forme=='super' else 'ellipse '}{r.param:g}"
g['nom']=g.apply(lab,axis=1)
# fenêtres
win=[(100,300),(300,1000),(1000,3000),(3000,10000),(10000,30001)]
print("C_eff par fenêtre (moyenne ± erreur, sur 64 échantillons × R de la fenêtre)")
hdr="forme".ljust(16)+"".join(f"[{a},{b-1 if b>30000 else b}]".rjust(17) for a,b in win)
print(hdr)
tab={}
for nom,s in d.assign(nom=d.apply(lab,axis=1)).groupby('nom',sort=False):
    row=[]
    for a,b in win:
        x=s[(s.R>=a)&(s.R<b)]
        # moyenne par R puis moyenne des R (erreur : sem des échantillons groupés)
        m=x.groupby('R').C.mean(); se=np.sqrt((x.groupby('R').C.sem()**2).sum())/len(m)
        row.append((m.mean(),se))
    tab[nom]=row
    print(nom.ljust(16)+"".join(f"{m:9.4f}±{se:.4f}".rjust(17) for m,se in row))
print("\nDispersion entre formes (CV) par fenêtre :")
for j,(a,b) in enumerate(win):
    v=np.array([tab[n][j][0] for n in tab]); se=np.mean([tab[n][j][1] for n in tab])
    print(f"  [{a},{b}] moyenne {v.mean():.4f}  CV {100*v.std()/v.mean():.2f} %   bruit attendu ~{100*se/v.mean():.2f} %")
print("\nNombre de sommets m/(Ω R^2/3), R >= 1000 :",
      round(d[d.R>=1000].cm.mean(),4), "±", round(d[d.R>=1000].groupby('nom' if 'nom' in d else ['forme','param']).cm.mean().std(),4))
# dérive globale
G=d.groupby('R').C.mean()
for x,lbl in [(np.log(G.index),'ln R'),(G.index**(-1/3),'R^-1/3')]:
    m=G.index>=300
    c=np.polyfit(x[m],G.values[m],1); print(f"ajustement linéaire C vs {lbl} (R≥300) : pente {c[0]:.4f}, ordonnée {c[1]:.4f}")
g.to_csv('C_eff_par_R.csv',index=False)
print("\nC moyen (11 formes) par R :"); print(G.round(4).to_string())
