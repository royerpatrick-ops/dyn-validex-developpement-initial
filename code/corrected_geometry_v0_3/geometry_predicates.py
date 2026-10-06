"""Closed-body lattice membership, with exact rationals or filtered intervals.

Inputs int/Fraction are exact. A float means its exact binary rational value,
including theta (whose mathematical sin/cos, not libm-rounded values, are used).
No membership tolerance is used. Ambiguity at the precision limit raises.
The interval path assumes IEEE round-to-nearest basic NumPy operations and
correct mpmath 1.4.1 interval elementary functions; it is subject to audit.
"""
from fractions import Fraction
import math
import numpy as np
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import to_rational

class GeometryIndeterminate(ArithmeticError):
    pass

def rational(value):
    if isinstance(value, bool):
        raise ValueError('boolean is not a geometry parameter')
    if isinstance(value, Fraction):
        result = value
    elif isinstance(value, (int, np.integer)):
        result = Fraction(int(value))
    elif isinstance(value, (float, np.floating)) and math.isfinite(float(value)):
        result = Fraction(float(value))
    else:
        raise ValueError('finite int, float or Fraction required')
    if not math.isfinite(float(result)):
        raise ValueError('finite parameter required')
    return result

def endpoints(value):
    return tuple(Fraction(*to_rational(endpoint)) for endpoint in value._mpi_)

def iv_rational(ctx, value):
    return ctx.mpf(value.numerator) / ctx.mpf(value.denominator)

def outward_rational(value):
    f = float(value)
    represented = Fraction(f)
    return (math.nextafter(f, -math.inf) if represented > value else f,
            math.nextafter(f, math.inf) if represented < value else f)

def integer_arrays(kx,ky):
    kx,ky=np.broadcast_arrays(np.asarray(kx),np.asarray(ky))
    if not all(np.issubdtype(a.dtype,np.integer) for a in (kx,ky)):
        raise ValueError('integer lattice indices required; floats are not truncated')
    if any(np.any((a < -2**40) | (a > 2**40)) for a in (kx,ky)):
        raise ValueError('indices exceed the exact float64-filter range')
    return kx.astype(np.int64),ky.astype(np.int64)

class Interval:
    """Float64 interval; every basic operation expands its rounded endpoints."""
    def __init__(self, lo, hi=None):
        self.lo = np.asarray(lo, dtype=float)
        self.hi = self.lo if hi is None else np.asarray(hi, dtype=float)

    @staticmethod
    def enclosing(lo, hi):
        bad = np.isnan(lo) | np.isnan(hi)
        lo = np.where(bad, -np.inf, np.nextafter(lo, -np.inf))
        hi = np.where(bad, np.inf, np.nextafter(hi, np.inf))
        return Interval(lo, hi)

    def __add__(self, other):
        with np.errstate(all='ignore'):
            return self.enclosing(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        with np.errstate(all='ignore'):
            products = np.array([self.lo*other.lo, self.lo*other.hi,
                                 self.hi*other.lo, self.hi*other.hi])
            return self.enclosing(np.min(products,axis=0),np.max(products,axis=0))

    def __truediv__(self, other):
        if np.any((other.lo <= 0) & (other.hi >= 0)):
            return Interval(np.full(np.broadcast_shapes(self.lo.shape,other.lo.shape),-np.inf),
                            np.full(np.broadcast_shapes(self.hi.shape,other.hi.shape),np.inf))
        with np.errstate(all='ignore'):
            inverse = self.enclosing(1/other.hi, 1/other.lo)
            return self * inverse

    def absolute(self):
        lo = np.where(self.lo>=0,self.lo,np.where(self.hi<=0,-self.hi,0))
        return Interval(lo,np.maximum(np.abs(self.lo),np.abs(self.hi)))

    def integer_power(self, n):
        # Caller uses a nonnegative interval; repeated products avoid np.power's
        # unspecified elementary-function error at nontrivial exponents.
        result = Interval(np.ones_like(self.lo))
        base = self
        while n:
            if n & 1:
                result = result * base
                result.lo = np.maximum(result.lo,0)
            n //= 2
            if n:
                base = base * base
                base.lo = np.maximum(base.lo,0)
        return result

    def positive_power(self, exponent):
        if exponent.denominator == 1:
            return self.integer_power(exponent.numerator)
        if exponent.denominator == 2:
            with np.errstate(all='ignore'):
                root = self.enclosing(np.sqrt(np.maximum(self.lo,0)),np.sqrt(np.maximum(self.hi,0)))
            root.lo = np.maximum(root.lo,0)
            return root.integer_power(exponent.numerator)
        # Non-half-integer powers bypass the float64 filter entirely.
        return Interval(np.full_like(self.lo,-np.inf),np.full_like(self.hi,np.inf))

def exact_power(base, exponent):
    if base in (0,1):
        return base
    if exponent.denominator == 1:
        return base ** exponent.numerator
    if exponent.denominator == 2:
        a,b = math.isqrt(base.numerator),math.isqrt(base.denominator)
        if a*a == base.numerator and b*b == base.denominator:
            return Fraction(a,b) ** exponent.numerator
    return None

class BodyPredicate:
    def __init__(self, kind, prm, R, theta, u=0, v=0, *, precisions=(50,100,200,400)):
        if kind not in ('ellipse','super'):
            raise ValueError('unknown shape')
        self.kind=kind
        self.prm,self.R,self.theta,self.u,self.v=map(rational,(prm,R,theta,u,v))
        if self.R<=0 or self.prm<=0 or (kind=='super' and self.prm<=1):
            raise ValueError('R>0, ellipse ratio>0 or superellipse p>1 required')
        if not (0<=self.u<1 and 0<=self.v<1):
            raise ValueError('translations must be in [0,1)')
        if not precisions or any(not isinstance(d,int) or d<15 for d in precisions):
            raise ValueError('positive precision schedule required')
        self.precisions=tuple(precisions)
        self.contexts={}
        self.stats={'filtered_decisions':0,'exact_decisions':0,'interval_refinements':0,
                    'indeterminate':0,'exact_symmetry_ties':0}
        self.disk=(kind=='ellipse' and self.prm==1) or (kind=='super' and self.prm==2)
        self.Rf=Interval(*outward_rational(self.R))
        self.uf=Interval(*outward_rational(self.u)); self.vf=Interval(*outward_rational(self.v))
        self.pf=Interval(*outward_rational(self.prm))
        if self.theta==0 or self.disk:
            self.cf,self.sf=Interval(1),Interval(0)
        else:
            ctx,c,s=self._context(self.precisions[0])
            cl,ch=endpoints(c);sl,sh=endpoints(s)
            self.cf=Interval(outward_rational(cl)[0],outward_rational(ch)[1])
            self.sf=Interval(outward_rational(sl)[0],outward_rational(sh)[1])

    def _context(self,dps):
        if dps not in self.contexts:
            ctx=MPIntervalContext();ctx.dps=dps
            t=iv_rational(ctx,self.theta)
            self.contexts[dps]=(ctx,ctx.cos(t),ctx.sin(t))
        return self.contexts[dps]

    @property
    def enclosing_radius(self):
        # Rotation preserves Euclidean norm. Ellipse lies in radius max(R,R/k);
        # any superellipse here lies in [-R,R]^2 and hence in radius 2R.
        return max(self.R,self.R/self.prm) if self.kind=='ellipse' else 2*self.R

    def bounds(self,kx,ky):
        x=(Interval(kx)-self.uf)/self.Rf
        y=(Interval(ky)-self.vf)/self.Rf
        X=(self.cf*x+self.sf*y).absolute()
        Y=((-self.sf)*x+self.cf*y).absolute()
        if self.kind=='ellipse':
            return X.integer_power(2)+(self.pf*Y).absolute().integer_power(2)
        return X.positive_power(self.prm)+Y.positive_power(self.prm)

    def _exact_value(self,x,y):
        if self.disk:
            return (x*x+y*y)/(self.R*self.R)
        if self.theta!=0:
            return None
        X,Y=abs(x/self.R),abs(y/self.R)
        if self.kind=='ellipse':
            return X*X+(self.prm*Y)**2
        a,b=exact_power(X,self.prm),exact_power(Y,self.prm)
        return None if a is None or b is None else a+b

    def _iv_value(self,x,y,dps):
        ctx,c,s=self._context(dps)
        xx,yy=iv_rational(ctx,x)/iv_rational(ctx,self.R),iv_rational(ctx,y)/iv_rational(ctx,self.R)
        X,Y=abs(c*xx+s*yy),abs(-s*xx+c*yy)
        if self.kind=='ellipse':
            return X**2+(iv_rational(ctx,self.prm)*Y)**2
        p=iv_rational(ctx,self.prm)
        return X**p+Y**p

    def contains_point(self,x,y):
        x,y=rational(x),rational(y)
        value=self._exact_value(x,y)
        if value is not None:
            self.stats['exact_decisions']+=1
            return value<=1
        for dps in self.precisions:
            self.stats['interval_refinements']+=1
            lo,hi=endpoints(self._iv_value(x,y,dps))
            if hi<=1:return True
            if lo>1:return False
        self.stats['indeterminate']+=1
        raise GeometryIndeterminate(f'membership unresolved at {self.precisions[-1]} digits: x={x}, y={y}')

    def contains_indices(self,kx,ky):
        kx,ky=integer_arrays(kx,ky)
        shape=kx.shape;kx=kx.ravel();ky=ky.ravel()
        result=np.empty(kx.size,dtype=bool)
        value=self.bounds(kx,ky)
        inside=value.hi<=1; outside=value.lo>1
        result[inside]=True;result[outside]=False
        self.stats['filtered_decisions']+=int(np.count_nonzero(inside|outside))
        for i in np.flatnonzero(~(inside|outside)):
            result[i]=self.contains_point(Fraction(int(kx[i]))-self.u,Fraction(int(ky[i]))-self.v)
        return result.reshape(shape)

    def _compare_exact_or_iv(self,kx,ky0,ky1):
        x=Fraction(int(kx))-self.u
        y0,y1=Fraction(int(ky0))-self.v,Fraction(int(ky1))-self.v
        if y0==-y1 and (x==0 or self.theta==0):
            self.stats['exact_symmetry_ties']+=1
            return 0
        a,b=self._exact_value(x,y0),self._exact_value(x,y1)
        if a is not None and b is not None:
            self.stats['exact_decisions']+=1
            return (a>b)-(a<b)
        for dps in self.precisions:
            self.stats['interval_refinements']+=1
            lo,hi=endpoints(self._iv_value(x,y0,dps)-self._iv_value(x,y1,dps))
            if hi<0:return -1
            if lo>0:return 1
            if lo==hi==0:return 0
        self.stats['indeterminate']+=1
        raise GeometryIndeterminate('discrete convex comparison unresolved; no column is guessed')

    def strictly_decreasing_step(self,kx,ky):
        # True iff f(ky)>f(ky+1): discrete minimum lies to the right.
        # Equality chooses the leftmost minimum.
        kx,ky=integer_arrays(kx,ky)
        difference=self.bounds(kx,ky)-self.bounds(kx,ky+1)
        decreasing=difference.lo>0; increasing=difference.hi<0
        result=decreasing.copy()
        self.stats['filtered_decisions']+=int(np.count_nonzero(decreasing|increasing))
        for i in np.flatnonzero(~(decreasing|increasing)):
            result.ravel()[i]=self._compare_exact_or_iv(int(kx.ravel()[i]),int(ky.ravel()[i]),int(ky.ravel()[i])+1)>0
        return result
