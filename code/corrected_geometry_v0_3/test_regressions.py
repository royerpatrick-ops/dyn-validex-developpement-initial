"""Targeted regressions for corrected mathematical/numerical failure modes."""
import ast
import csv
import hashlib
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch
import warnings
import numpy as np
from scipy.integrate import IntegrationWarning
from scipy.spatial import ConvexHull
import invariant_checks as inv
import integral_model as model
import formule_integrale as integral
import test_classes_w0 as classes
import balayage_translations as bt
from coefficient_B import cutoff_integral, model_coefficients

ROOT = Path(__file__).resolve().parent.parent

def integer_hull(points):
    points = sorted(set(map(tuple,points)))
    def cross(a,b,c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def half(seq):
        h=[]
        for p in seq:
            while len(h)>=2 and cross(h[-2],h[-1],p)<=0: h.pop()
            h.append(p)
        return h[:-1]
    return set(half(points)+half(points[::-1]))

class CorrectiveRegression(unittest.TestCase):
    def test_divergence_is_rejected_before_quadrature(self):
        with patch.object(inv,'quad',side_effect=AssertionError('quad must not run')):
            for p,a in [(1.3,1.5),(1.25,4/3),(4/3,1.5),(1.2,4/3)]:
                with self.subTest(p=p,a=a), self.assertRaises(inv.DivergentInvariant):
                    inv.I_super(p,a)

    def test_finite_counterexample_I43(self):
        self.assertAlmostEqual(inv.I_super(1.3,4/3),22.138237169,places=7)

    def test_disk_invariants(self):
        for a in [0,1/3,4/3,1.5]:
            self.assertAlmostEqual(inv.I_super(2,a),2*math.pi,places=11)
            self.assertAlmostEqual(inv.I_ellipse(1,a),2*math.pi,places=11)

    def test_supplied_invariants_preserved(self):
        with (ROOT/'originals/invariants.csv').open(newline='') as f:
            for r in csv.DictReader(f):
                values = inv.invariants(r['forme'],float(r['param']))
                np.testing.assert_allclose(values,[float(r[k]) for k in ['P','Omega','I43']],rtol=1e-10,atol=1e-10)

    def test_warning_is_fatal(self):
        def bad(*args,**kwargs):
            warnings.warn('probably divergent',IntegrationWarning)
            return -28.543699117,1.
        with patch.object(inv,'quad',side_effect=bad), self.assertRaises(IntegrationWarning):
            inv.I_ellipse(2,1.5)

    def test_invalid_quadrature_is_fatal(self):
        for value,error in [(-1,0),(float('nan'),0),(1,float('inf')),(1,-1),(1,.1)]:
            with patch.object(inv,'quad',return_value=(value,error)), self.assertRaises(ArithmeticError):
                inv.I_ellipse(2,1.5)

    def test_invalid_invariant_parameters(self):
        for args in [('wrong',2,1),('super',1,1),('ellipse',0,1),('super',float('nan'),1),('ellipse',2,-1)]:
            with self.assertRaises(ValueError): inv.parameters(*args)

    def test_correct_class_coefficient(self):
        self.assertAlmostEqual(classes.K,4*math.sqrt(2)/15,places=15)
        self.assertEqual(classes.primitive_counts(5),{1:4,2:4,5:8})
        pred=4*classes.K*100000**(-1/6)/(2*math.pi)
        self.assertAlmostEqual(pred,.03523956345289601,places=14)

    def test_finite_tail_with_nonconstant_a3G(self):
        A,M,L,b = 6.,16.,1/9,3.
        exact=L*(1/A-1/M)+b/3*(A**-3-M**-3)
        value=integral.deterministic_tail(lambda a:L/a**3+b/a**5,A,M)
        self.assertAlmostEqual(value,exact,places=14)
        wrong=(L+b/((A+M)/2)**2)*(1/A-1/M)
        self.assertGreater(abs(wrong-exact),1e-4)

    def test_sampler_cutoffs_and_inputs(self):
        for a,n,kmax in [(0,2,1),(1,1,1),(float('nan'),2,1),(1,2,0)]:
            with self.assertRaises(ValueError): model.echantillon(a,n,np.random.default_rng(7),kmax)
        with self.assertRaises(RuntimeError): model.echantillon(10,100,np.random.default_rng(7),1)
        with self.assertRaises(RuntimeError): model.G_cible(10,100,np.random.default_rng(7),kmax=1)
        with self.assertRaises(ValueError): model.G_cible(1,100,np.random.default_rng(7))
        with self.assertRaises(ValueError): integral.estimate(A=6,Amax=5)

    def test_sampler_independent_scalar(self):
        for a in [.1,1.,4.,8.]:
            seed,n=37,200
            rng=np.random.default_rng(seed)
            delta=rng.uniform(0,1/a,n);s=rng.uniform(0,a,n);beta=rng.uniform(0,a,n)
            c,flags=model.echantillon(a,n,np.random.default_rng(seed))
            expected=[];edges=[]
            for d0,s0,b0 in zip(delta,s,beta):
                for k in range(20000):
                    d=d0+k/a; off=(s0+k*b0)%a;r=math.sqrt(2*d)
                    js=[j for j in range(math.floor((-r-off)/a)-1,math.ceil((r-off)/a)+2) if abs(off+j*a)<=r]
                    if js:
                        edges.append(len(js)>=2)
                        if len(js)>=2:
                            x1,x2=off+min(js)*a,off+max(js)*a
                            expected.append(d*(x2-x1)-(x2**3-x1**3)/3)
                        else: expected.append(0.)
                        break
                else: self.fail('scalar sampler unresolved')
            np.testing.assert_allclose(c,expected,atol=1e-12,rtol=1e-12)
            np.testing.assert_array_equal(flags,edges)

    def test_exact_rational_disk_vertices(self):
        for R in [5,30,1000]:
            for numerator_u,numerator_v,den in [(0,0,1),(1,1,2),(1,2,4),(3,7,8)]:
                exact=[]
                for kx in range(-R-1,R+2):
                    xnum=den*kx-numerator_u
                    square=(den*R)**2-xnum*xnum
                    if square<0: continue
                    h=math.isqrt(square)
                    kb=-((-numerator_v+h)//den)
                    kt=(numerator_v+h)//den
                    if kb<=kt: exact.extend([(kx,kb),(kx,kt)])
                pts=bt.columns('ellipse',1.,float(R),.31,numerator_u/den,numerator_v/den)
                got=np.rint(pts[ConvexHull(pts).vertices]+[numerator_u/den,numerator_v/den]).astype(int)
                self.assertEqual(set(map(tuple,got)),integer_hull(exact))

    def test_perimeter_requires_integer_indices(self):
        square = np.array([[0,0],[1,0],[1,1],[0,1]],dtype=np.int64)
        self.assertEqual(bt.hull_perim(square),(4.,4))
        with self.assertRaises(ValueError): bt.hull_perim(square.astype(float))

    def test_B_finite_cutoff_independent_sum(self):
        T=8; lam=1/T; contributions=[]
        for x in range(-T,T+1):
            for y in range(-T,T+1):
                r=math.hypot(x,y)
                if r and r<T and math.gcd(x,y)==1:
                    contributions.append((8/5)*r**-1.5-2*math.sqrt(lam)/r+(2/5)*lam**2.5*r)
        expected=math.pi/8*math.fsum(contributions)-2/math.sqrt(lam)
        self.assertAlmostEqual(cutoff_integral(T),expected,places=12)
        self.assertAlmostEqual(model_coefficients()['B_total_model'],.2894234178601282,places=13)

    def test_original_sources_intact(self):
        manifest=json.loads((ROOT/'originals_sha256.json').read_text())
        self.assertEqual(len(manifest),89)
        actual={str(p.relative_to(ROOT/'originals')) for p in (ROOT/'originals').rglob('*') if p.is_file()}
        self.assertEqual(actual,set(manifest))
        for name,h in manifest.items():
            self.assertEqual(hashlib.sha256((ROOT/'originals'/name).read_bytes()).hexdigest(),h)

if __name__ == '__main__':
    unittest.main(verbosity=2)
