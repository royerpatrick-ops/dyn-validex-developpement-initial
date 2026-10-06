"""Column extrema from exact rational bounds and decisive membership predicates.

No continuous minimizer, discriminant or floating floor/ceil determines exclusion.
The non-disk search uses strict convexity to find a discrete minimum, then two
monotone searches. Ambiguous predicates abort the call. 'certified' refers to
this algorithm's acceptance rule, not a certification of the asymptotic law.
"""
from fractions import Fraction
import math
import numpy as np
from geometry_predicates import BodyPredicate, rational

def floor_fraction(x):
    return x.numerator//x.denominator

def ceil_fraction(x):
    return -((-x.numerator)//x.denominator)

def disk_column_indices(R,u=0,v=0,*,max_columns=2000005):
    R,u,v=map(rational,(R,u,v))
    if R<=0 or not (0<=u<1 and 0<=v<1):
        raise ValueError('R>0 and translations in [0,1) required')
    xmin,xmax=ceil_fraction(u-R),floor_fraction(u+R)
    if xmax-xmin+1>max_columns:
        raise ValueError('column budget exceeded')
    D=math.lcm(R.denominator,u.denominator,v.denominator)
    radius=int(R*D);U=int(u*D);V=int(v*D)
    points=[]
    for kx in range(xmin,xmax+1):
        square=radius*radius-(D*kx-U)**2
        if square<0:continue
        h=math.isqrt(square)
        bottom=-((h-V)//D);top=(h+V)//D
        if bottom<=top:
            points.extend(((kx,bottom),(kx,top)))
    if points and max(abs(z) for pair in points for z in pair)>2**40:
        raise ValueError('indices exceed supported range')
    return np.asarray(points,dtype=np.int64).reshape(-1,2)

def _generic_batch(body,kx,ylo,yhi):
    lo=np.full(kx.size,ylo,dtype=np.int64);hi=np.full(kx.size,yhi,dtype=np.int64)
    while np.any(lo<hi):
        active=np.flatnonzero(lo<hi);mid=(lo[active]+hi[active])//2
        decrease=body.strictly_decreasing_step(kx[active],mid)
        lo[active[decrease]]=mid[decrease]+1
        hi[active[~decrease]]=mid[~decrease]
    minimum=lo
    good=body.contains_indices(kx,minimum)
    kx,minimum=kx[good],minimum[good]
    if not kx.size:return np.empty((0,2),dtype=np.int64)
    lo=np.full(kx.size,ylo,dtype=np.int64);hi=minimum.copy()
    while np.any(lo<hi):
        active=np.flatnonzero(lo<hi);mid=(lo[active]+hi[active])//2
        inside=body.contains_indices(kx[active],mid)
        hi[active[inside]]=mid[inside]
        lo[active[~inside]]=mid[~inside]+1
    bottom=lo.copy()
    lo=minimum.copy();hi=np.full(kx.size,yhi,dtype=np.int64)
    while np.any(lo<hi):
        active=np.flatnonzero(lo<hi);mid=(lo[active]+hi[active]+1)//2
        inside=body.contains_indices(kx[active],mid)
        lo[active[inside]]=mid[inside]
        hi[active[~inside]]=mid[~inside]-1
    top=lo
    # Assert endpoint and neighbour decisions: guards the search implementation.
    if not (body.contains_indices(kx,bottom).all() and body.contains_indices(kx,top).all()
            and not body.contains_indices(kx,bottom-1).any() and not body.contains_indices(kx,top+1).any()):
        raise ArithmeticError('column extrema verification failed')
    return np.vstack((np.column_stack((kx,bottom)),np.column_stack((kx,top))))

def column_indices(kind,prm,R,theta,u=0,v=0,*,max_columns=2000005,batch_size=4096,precisions=(50,100,200,400),return_stats=False):
    if not isinstance(max_columns,int) or max_columns<1 or not isinstance(batch_size,int) or batch_size<1:
        raise ValueError('positive integer budgets required')
    body=BodyPredicate(kind,prm,R,theta,u,v,precisions=precisions)
    if body.disk:
        result=disk_column_indices(body.R,body.u,body.v,max_columns=max_columns)
        stats={'path':'EXACT_RATIONAL_DISK','returned_extrema':len(result)}
    else:
        B=body.enclosing_radius
        xmin,xmax=ceil_fraction(body.u-B),floor_fraction(body.u+B)
        ylo,yhi=floor_fraction(body.v-B),ceil_fraction(body.v+B)
        if max(abs(xmin),abs(xmax),abs(ylo),abs(yhi))>2**40-1:
            raise ValueError('bounding box exceeds supported integer range')
        if xmax-xmin+1>max_columns:
            raise ValueError('column budget exceeded')
        parts=[]
        for begin in range(xmin,xmax+1,batch_size):
            kx=np.arange(begin,min(begin+batch_size,xmax+1),dtype=np.int64)
            parts.append(_generic_batch(body,kx,ylo,yhi))
        result=np.vstack(parts) if parts else np.empty((0,2),dtype=np.int64)
        stats={'path':'FILTERED_CONVEX_INTEGER_SEARCH',**body.stats,'returned_extrema':len(result)}
    return (result,stats) if return_stats else result

def columns(kind,prm,R,theta,u=0,v=0):
    """Compatibility coordinates; membership is decided before conversion to float."""
    indices=column_indices(kind,prm,R,theta,u,v)
    return indices.astype(float)-[float(rational(u)),float(rational(v))]

def integer_hull_vertices(points):
    """Andrew hull with Python-integer orientation, no Qhull or float predicate."""
    array=np.asarray(points)
    if array.size==0:return np.empty((0,2),dtype=np.int64)
    if array.ndim!=2 or array.shape[1]!=2 or not np.issubdtype(array.dtype,np.integer):
        raise ValueError('an N by 2 array of integer indices is required')
    if np.any((array < -2**40) | (array > 2**40)):
        raise ValueError('integer coordinates outside supported range')
    points=sorted(set((int(x),int(y)) for x,y in array))
    if len(points)<3:
        return np.asarray(points,dtype=np.int64).reshape(-1,2)
    def cross(a,b,c):
        return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
    def chain(seq):
        result=[]
        for p in seq:
            while len(result)>=2 and cross(result[-2],result[-1],p)<=0:result.pop()
            result.append(p)
        return result[:-1]
    return np.asarray(chain(points)+chain(points[::-1]),dtype=np.int64).reshape(-1,2)

def integer_hull_perimeter(points):
    vertices=integer_hull_vertices(points)
    if len(vertices)<3:
        raise ValueError('nondegenerate hull required for the perimeter convention')
    pairs=list(zip(vertices,np.roll(vertices,-1,axis=0)))
    perimeter=math.fsum(math.hypot(int(b[0])-int(a[0]),int(b[1])-int(a[1])) for a,b in pairs)
    return perimeter,len(vertices)
