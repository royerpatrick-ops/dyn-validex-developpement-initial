"""Compact independent geometry evidence; no C calculation or new campaign.

Run from package root: python candidate/run_geometry_checks.py --out NEW_DIRECTORY
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import time
import numpy as np
import scipy
import mpmath
from geometry_columns import column_indices,integer_hull_perimeter
from test_geometry_v03 import reference_disk,vertex_set,ANGLES,SHIFTS,FAILED_RADII

ROOT=Path(__file__).resolve().parent.parent

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',required=True,type=Path)
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    clock=time.monotonic()
    result={'status':'PASS_LOCAL_GEOMETRY_REGRESSIONS',
            'scope':'Finite geometry regressions; no full numerical/asymptotic certification.',
            'C_recalculated':False,'full_campaign_run':False,
            'environment':{'python':platform.python_version(),'platform':platform.platform(),
                           'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mpmath.__version__},
            'candidate_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted((ROOT/'candidate').glob('*.py'))}}
    tested=0;failures=[]
    for R in range(2,251):
        expected=reference_disk(R)
        for angle in ANGLES:
            actual=vertex_set(column_indices('ellipse',1,R,angle))
            tested+=1
            if actual!=expected:failures.append({'R':R,'theta':angle,'missing':sorted(expected-actual),'extra':sorted(actual-expected)})
    result['all_radii_2_250']={'cases':tested,'failures':failures,'passed':tested-len(failures)}
    details=[]
    for R in [31,39,41,43,53,1000]:
        indices,stats=column_indices('ellipse',1,R,.314159,return_stats=True)
        expected=reference_disk(R)
        actual=vertex_set(indices)
        details.append({'R':R,'expected_vertices':len(expected),'actual_vertices':len(actual),
                        'missing':sorted(expected-actual),'extra':sorted(actual-expected),
                        'perimeter':integer_hull_perimeter(indices)[0],'predicate_path':stats['path']})
    result['previous_counterexamples']=details
    cases=0;wrong=[]
    for R in FAILED_RADII+[1000]:
        for u,v in SHIFTS:
            expected=reference_disk(R,u,v)
            for angle in ANGLES:
                cases+=1
                actual=vertex_set(column_indices('ellipse',1,R,angle,u,v))
                if actual!=expected:wrong.append({'R':R,'u':str(u),'v':str(v),'theta':angle})
    result['failed_radii_and_1000_with_shifts']={'cases':cases,'passed':cases-len(wrong),'failures':wrong}

    # Replay existing stored individual hulls, not a new Monte Carlo campaign.
    with (ROOT/'originals/balayage_translations.csv').open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    selected=[]
    for r in rows:
        kind,prm,R,i=r['forme'],float(r['param']),int(float(r['R'])),int(r['echantillon'])
        if i==0 and (R in [100,1230] or (R==30000 and (kind,prm) in [('ellipse',2.),('super',1.5),('super',10.)])):
            selected.append(r)
    replays=[]
    for r in selected:
        t=time.monotonic()
        indices,stats=column_indices(r['forme'],float(r['param']),float(r['R']),float(r['theta']),float(r['u']),float(r['v']),return_stats=True)
        perim,n=integer_hull_perimeter(indices)
        replays.append({'forme':r['forme'],'param':float(r['param']),'R':float(r['R']),
                        'sample':int(r['echantillon']),'perimeter_difference':perim-float(r['P_hull']),
                        'old_vertex_count':int(r['m_sommets']),'new_vertex_count':n,
                        'geometry_stats':stats,'seconds':time.monotonic()-t})
    result['old_csv_replays']=replays
    result['old_csv_replay_max_abs_perimeter_difference']=max(abs(r['perimeter_difference']) for r in replays)
    result['old_csv_replay_same_vertex_count']=all(r['old_vertex_count']==r['new_vertex_count'] for r in replays)
    result['seconds_total']=time.monotonic()-clock
    if failures or wrong or any(r['missing'] or r['extra'] for r in details):
        result['status']='FAIL_LOCAL_GEOMETRY_REGRESSIONS'
    (args.out/'result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    if result['status'].startswith('FAIL'):
        raise SystemExit(1)

if __name__=='__main__':
    main()
