"""Recompute all-form residuals, plus the combined four-form large-R subset.

Chi-square tail probabilities assume independent normal estimates with known
errors. They omit uncertainty in shared C,B and the finite-R remainder.
"""
import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from invariant_checks import I_super, I_ellipse

ROOT = Path(__file__).resolve().parent.parent

def load(name):
    with name.open(newline='', encoding='utf-8-sig') as f:
        rows = list(csv.DictReader(f))
    keys = set()
    groups = defaultdict(list)
    for row in rows:
        kind, prm, R, i = row['forme'], float(row['param']), float(row['R']), int(row['echantillon'])
        key = kind, prm, R, i
        if key in keys:
            raise ValueError(f'doublon dans {name.name}: {key}')
        keys.add(key)
        vals = {k:float(v) for k,v in row.items() if k != 'forme'}
        if not all(math.isfinite(v) for v in vals.values()):
            raise ValueError('donnee non finie')
        if abs(vals['P_cont']-vals['P_hull']-vals['deficit_abs']) > 1e-7:
            raise ValueError('identite de deficit rompue')
        groups[(kind,prm,R)].append(vals['deficit_abs'])
    return rows, groups

def summarize(groups):
    out = []
    for (kind, prm, R), vals in sorted(groups.items()):
        I = I_super if kind == 'super' else I_ellipse
        J, L = I(prm,4/3), I(prm,1.5)
        c = np.asarray(vals)*R**(1/3)/J
        se = float(c.std(ddof=1)/math.sqrt(len(c)))
        if not math.isfinite(se) or se <= 0:
            raise ValueError('erreur-type positive requise')
        prediction = .7192 - .2894*R**(-1/6)*L/J
        out.append({'forme':kind,'param':prm,'R':R,'count':len(c), 'mean':float(c.mean()),
                    'mc_se':se,'prediction':prediction,'z':(float(c.mean())-prediction)/se})
    return out

def score(rows):
    s = math.fsum(r['z']**2 for r in rows)
    return {'groups':len(rows), 'chi2':s, 'chi2_per_group':s/len(rows),
            'nominal_p_independent_normal_known_errors':float(chi2.sf(s,len(rows)))}

def run(source):
    r1, g1 = load(source/'balayage_translations.csv')
    r2, g2 = load(source/'balayage_grands_R.csv')
    allrows = summarize(g1)
    combined = defaultdict(list)
    shapes = {(k,p) for k,p,R in g2}
    for groups in [g1,g2]:
        for (k,p,R), vals in groups.items():
            if (k,p) == ('super',2.): k,p = 'ellipse',1.
            if (k,p) in shapes: combined[(k,p,R)].extend(vals)
    four = [r for r in summarize(combined) if r['R']>=1230]
    return {'status':'DESCRIPTIVE_FINITE_R_DIAGNOSTIC',
            'constants_fixed_for_comparison':{'C':.7192,'B':.2894},
            'assumptions':'Independent normal means with known SEM; shared constants and finite-R remainder omitted.',
            'source_counts':{ 'campaign_1_rows':len(r1),'campaign_2_rows':len(r2)},
            'all_11_R_ge_1230':score([r for r in allrows if r['R']>=1230]),
            'all_11_R_ge_10000':score([r for r in allrows if r['R']>=10000]),
            'super_1_5_R_ge_1230':score([r for r in allrows if r['R']>=1230 and (r['forme'],r['param'])==('super',1.5)]),
            'four_forms_R_ge_1230':score(four),
            'counts_at_R_30000':[r for r in four if r['R']==30000]}, allrows, four

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=ROOT/'originals')
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    result, allrows, four = run(args.source)
    (args.out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    for name, rows in [('all_forms.csv',allrows),('four_forms.csv',four)]:
        with (args.out/name).open('w',newline='',encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)
    print(json.dumps(result,indent=2))

if __name__ == '__main__':
    main()
