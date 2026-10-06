"""Geometry regressions: exact boundary, rational shifts, filtered predicates.

The brute reference uses full-grid enumeration and Jarvis wrapping, independently
of the production column search and Andrew hull. Large disk references use an
alternative Fraction residual per column, then Qhull on the integer endpoints.
"""
from fractions import Fraction as F
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from scipy.spatial import ConvexHull
from mpmath.ctx_mp import MPContext
from geometry_predicates import BodyPredicate, GeometryIndeterminate, endpoints, Interval
from geometry_columns import column_indices, disk_column_indices, integer_hull_vertices, integer_hull_perimeter

ROOT=Path(__file__).resolve().parent.parent
FAILED_RADII=[31,39,41,43,53,55,62,63,78,82,85,86,106,110,124,126,156,164,170,172,207,212,220,239,248]
ANGLES=[0.,.31,.314159,math.pi/4,math.pi/2,-.314159]
SHIFTS=[(F(0),F(0)),(F(1,2),F(1,2)),(F(1,4),F(1,2)),(F(3,8),F(7,8)),(F(1,3),F(2,3))]

def vertex_set(points):
    return {tuple(map(int,p)) for p in integer_hull_vertices(points)}

def jarvis(points):
    """Reference differs from production's monotone chain."""
    points=sorted(set(points))
    if len(points)<3:return set(points)
    start=points[0];current=start;vertices=[]
    for _ in range(len(points)+1):
        vertices.append(current)
        candidate=next(p for p in points if p!=current)
        for p in points:
            if p==current:continue
            cross=(candidate[0]-current[0])*(p[1]-current[1])-(candidate[1]-current[1])*(p[0]-current[0])
            distance=lambda a:(a[0]-current[0])**2+(a[1]-current[1])**2
            if cross<0 or (cross==0 and distance(p)>distance(candidate)):
                candidate=p
        current=candidate
        if current==start:return set(vertices)
    raise AssertionError('reference hull failed to close')

def reference_disk(R,u=F(0),v=F(0)):
    R=F(R);u=F(u);v=F(v)
    pts=[]
    for kx in range(math.ceil(u-R),math.floor(u+R)+1):
        residual=R*R-(F(kx)-u)**2
        if residual<0:continue
        # Rational sqrt bracket: isqrt(n*d)/d. Its width is <1/d and
        # candidates are checked against the original squared inequality.
        h=math.isqrt(residual.numerator*residual.denominator)
        root_lo=F(h,residual.denominator);root_hi=F(h+1,residual.denominator)
        bottom=math.floor(v-root_hi)
        top=math.ceil(v+root_hi)
        while (F(kx)-u)**2+(F(bottom)-v)**2>R*R:bottom+=1
        while (F(kx)-u)**2+(F(top)-v)**2>R*R:top-=1
        if bottom<=top:pts.extend([(kx,bottom),(kx,top)])
    array=np.asarray(pts,dtype=np.int64).reshape(-1,2)
    return {tuple(map(int,p)) for p in array[ConvexHull(array).vertices]}

def independent_full_grid(kind,prm,R,theta,u,v):
    ctx=MPContext();ctx.dps=120
    def number(q):
        q=F(q);return ctx.mpf(q.numerator)/q.denominator
    radius,parameter,t=map(number,(R,prm,theta))
    c,s=ctx.cos(t),ctx.sin(t)
    B=math.ceil(max(F(R),F(R)/F(prm))) if kind=='ellipse' else math.ceil(2*F(R))
    points=[]
    for kx in range(-B-1,B+2):
        for ky in range(-B-1,B+2):
            x,y=F(kx)-F(u),F(ky)-F(v)
            if (kind=='ellipse' and F(prm)==1) or (kind=='super' and F(prm)==2):
                inside=x*x+y*y<=F(R)**2
            elif F(theta)==0 and (kind=='ellipse' or F(prm).denominator==1):
                X,Y=abs(x/F(R)),abs(y/F(R))
                value=X*X+(F(prm)*Y)**2 if kind=='ellipse' else X**int(prm)+Y**int(prm)
                inside=value<=1
            else:
                X,Y=abs((c*number(x)+s*number(y))/radius),abs((-s*number(x)+c*number(y))/radius)
                value=X**2+(parameter*Y)**2 if kind=='ellipse' else X**parameter+Y**parameter
                # Axis endpoints for noninteger p are exactly known, not rounded.
                if F(theta)==0 and (abs(x)==F(R) and y==0 or abs(y)==F(R) and x==0):
                    inside=True
                else:
                    if abs(value-1)<ctx.mpf('1e-100'):
                        raise AssertionError('independent reference ambiguous')
                    inside=value<1
            if inside:points.append((kx,ky))
    return jarvis(points)

class GeometryV03(unittest.TestCase):
    def test_all_249_radii_multiple_angles(self):
        for R in range(2,251):
            expected=reference_disk(R)
            for th in ANGLES:
                with self.subTest(R=R,theta=th):
                    self.assertEqual(vertex_set(column_indices('ellipse',1,R,th)),expected)

    def test_25_failures_rational_shifts_angles(self):
        for R in FAILED_RADII+[1000]:
            for u,v in SHIFTS:
                expected=reference_disk(R,u,v)
                for th in ANGLES:
                    with self.subTest(R=R,u=u,v=v,theta=th):
                        self.assertEqual(vertex_set(column_indices('ellipse',1,R,th,u,v)),expected)

    def test_original_24_extended_cases(self):
        for R in [5,30,1000]:
            for u,v in SHIFTS[:4]:
                for th in [.31,.314159]:
                    self.assertEqual(vertex_set(column_indices('ellipse',1,R,th,u,v)),reference_disk(R,u,v))

    def test_disk_exhaustive_full_grid_independent(self):
        for R in [F(7,2),F(5),F(11,3)]:
            for u,v in SHIFTS:
                points=[]
                for x in range(-7,8):
                    for y in range(-7,8):
                        if (F(x)-u)**2+(F(y)-v)**2<=R*R:points.append((x,y))
                self.assertEqual(vertex_set(disk_column_indices(R,u,v)),jarvis(points))

    def test_boundary_neighbours_exact(self):
        for R,x,y in [(31,31,0),(39,36,15),(41,40,9),(1000,800,600)]:
            predicate=BodyPredicate('ellipse',1,R,.314159)
            self.assertTrue(predicate.contains_point(x,y))
            self.assertTrue(predicate.contains_point(math.nextafter(float(x),0),y))
            self.assertFalse(predicate.contains_point(math.nextafter(float(x),math.inf),y))
        # Exact rational translations and fractional radii are retained as exact.
        p=BodyPredicate('ellipse',1,F(31,3),.314159,F(1,3),F(2,3))
        self.assertTrue(p.contains_point(F(31,3),0))
        self.assertFalse(p.contains_point(F(31,3)+F(1,10**80),0))

    def test_radius_ulp_and_subnormal_shift(self):
        for R in [math.nextafter(31.,0),31.,math.nextafter(31.,math.inf)]:
            for u,v in [(0,0),(np.nextafter(0.,1.),0)]:
                self.assertEqual(vertex_set(column_indices('ellipse',1,R,.314159,u,v)),reference_disk(F(R),F(float(u)),F(v)))

    def test_superellipse_p2_is_exact_disk(self):
        self.assertEqual(vertex_set(column_indices('super',2,31,.314159,F(1,3),F(2,3))),reference_disk(31,F(1,3),F(2,3)))

    def test_all_campaign_shape_families_brute_grid(self):
        cases=[('super',p) for p in [1.5,2.5,3,4,5,10]]+[('ellipse',k) for k in [1,1.5,2,3,5]]
        for kind,prm in cases:
            for th,u,v in [(0,F(0),F(0)),(.314159,F(0),F(0)),(.314159,F(1,3),F(2,3))]:
                with self.subTest(kind=kind,prm=prm,theta=th):
                    expected=independent_full_grid(kind,prm,7,th,u,v)
                    self.assertEqual(vertex_set(column_indices(kind,prm,7,th,u,v)),expected)

    def test_non_half_integer_power_fallback(self):
        case=('super',F(13,10),5,.314159,F(1,4),F(3,4))
        indices,stats=column_indices(*case,return_stats=True)
        self.assertGreater(stats['interval_refinements'],0)
        self.assertEqual(stats['filtered_decisions'],0)
        self.assertEqual(vertex_set(indices),independent_full_grid(*case))

    def test_exact_axis_ellipse_with_noninteger_ratio(self):
        case=('ellipse',F(3,2),15,0,F(0),F(0))
        indices,stats=column_indices(*case,return_stats=True)
        self.assertGreater(stats['exact_decisions'],0)
        self.assertEqual(vertex_set(indices),independent_full_grid(*case))

    def test_symmetric_half_shift_ties(self):
        _,stats=column_indices('super',F(3,2),9,0,F(1,2),F(1,2),return_stats=True)
        self.assertGreater(stats['exact_symmetry_ties'],0)

    def test_float_filter_encloses_high_precision_values(self):
        rng=np.random.default_rng(42030)
        for kind,prm in [('ellipse',F(3,2)),('super',F(3,2)),('super',F(5,2)),('super',10)]:
            body=BodyPredicate(kind,prm,F(31,3),F(float(.314159)),F(1,3),F(2,3))
            x,y=rng.integers(-20,21,(2,80))
            bound=body.bounds(x,y)
            for i in range(len(x)):
                lo,hi=endpoints(body._iv_value(F(int(x[i]))-body.u,F(int(y[i]))-body.v,150))
                self.assertLessEqual(F(float(bound.lo[i])),lo)
                self.assertGreaterEqual(F(float(bound.hi[i])),hi)

    def test_indeterminate_does_not_accept_or_reject_silently(self):
        ctx=MPContext();ctx.dps=100;t=ctx.mpf(1)/10
        x,y=F(str(ctx.cos(t))),F(str(ctx.sin(t)))
        body=BodyPredicate('ellipse',2,1,F(1,10),precisions=(15,))
        with self.assertRaises(GeometryIndeterminate):body.contains_point(x,y)
        with patch.object(BodyPredicate,'strictly_decreasing_step',side_effect=GeometryIndeterminate('forced')):
            with self.assertRaises(GeometryIndeterminate):column_indices('ellipse',2,9,.314159)

    def test_resource_budget_refuses_before_partial_return(self):
        for kind,prm in [('ellipse',1),('ellipse',2),('super',3)]:
            with self.assertRaises(ValueError):column_indices(kind,prm,100,.314159,max_columns=2)
        with self.assertRaises(ValueError):column_indices('wrong',1,10,0)
        with self.assertRaises(ValueError):column_indices('ellipse',1,10,float('nan'))
        with self.assertRaises(ValueError):column_indices('ellipse',1,10,0,u=1)
        with self.assertRaises(ValueError):BodyPredicate('ellipse',2,10,0).contains_indices([1.2],[0])
        with self.assertRaises(ValueError):integer_hull_vertices([[0.1,0.],[1.,0.],[0.,1.]])

    def test_nonempty_and_empty_columns_against_full_grid(self):
        # Independent exact polynomials at theta=0; compares all column extrema,
        # so an error hidden strictly inside the final convex hull is detected.
        for kind,prm in [('ellipse',F(3,2)),('super',3),('super',4)]:
            R,u,v=F(9,2),F(1,3),F(2,3)
            expected={}
            for x in range(-10,11):
                ys=[]
                for y in range(-10,11):
                    X,Y=abs((F(x)-u)/R),abs((F(y)-v)/R)
                    value=X*X+(prm*Y)**2 if kind=='ellipse' else X**prm+Y**prm
                    if value<=1:ys.append(y)
                if ys:expected[x]=(min(ys),max(ys))
            actual={}
            for x,y in column_indices(kind,prm,R,0,u,v):
                x,y=int(x),int(y)
                low,high=actual.get(x,(y,y));actual[x]=(min(low,y),max(high,y))
            self.assertEqual(actual,expected)

    def test_integer_hull_does_not_overflow_orientation(self):
        R=2**40
        points=np.array([[-R,-R],[R,-R],[R,R],[-R,R],[0,0]],dtype=np.int64)
        self.assertEqual(vertex_set(points),{(-R,-R),(R,-R),(R,R),(-R,R)})
        self.assertEqual(integer_hull_perimeter(points),(float(8*R),4))
        with self.assertRaises(ValueError):integer_hull_perimeter(np.array([[0,0],[1,1]],dtype=np.int64))

    def test_subnormal_arithmetic_not_flushed(self):
        self.assertGreater(float(np.nextafter(0.,1.)),0)
        self.assertGreater(float(np.sqrt(np.nextafter(0.,1.))),0)
        self.assertGreater(float(np.float64(1e-160)*np.float64(1e-160)),0)

    def test_numeric_campaign_and_C_cli_are_suspended(self):
        with tempfile.TemporaryDirectory() as tmp:
            for script,args in [('balayage_translations.py',['--campaign','--out',str(Path(tmp)/'campaign')]),
                                ('formule_integrale.py',['--profile','reference','--out',str(Path(tmp)/'C')]),
                                ('test_classes_w0.py',['--out',str(Path(tmp)/'classes.json')])]:
                result=subprocess.run([sys.executable,str(ROOT/'candidate'/script),*args],capture_output=True,text=True)
                self.assertNotEqual(result.returncode,0)
                self.assertIn('suspendu',result.stderr)
                self.assertEqual(list(Path(tmp).iterdir()),[])

if __name__=='__main__':
    unittest.main(verbosity=2)
