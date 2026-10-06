import numpy as np, warnings, math
warnings.filterwarnings("ignore")
from scipy.integrate import quad
from scipy.spatial import ConvexHull
from balayage_corrige import hull_points, hull_perim, P_cont
from test_ellipses import pts_ellipse, invariants

def disk_pts(R, u=0.0, v=0.0):
    # disque de rayon R centré en (u,v) : colonnes exactes
    xs=np.arange(math.floor(u-R)-1, math.ceil(u+R)+2, dtype=float)
    d=R*R-(xs-u)**2; ok=d>=0; xs,d=xs[ok],d[ok]
    yt=np.floor(v+np.sqrt(d)+1e-12); yb=np.ceil(v-np.sqrt(d)-1e-12)
    f=lambda x,y:(x-u)**2+(y-v)**2
    yt[f(xs,yt)>R*R]-=1; yb[f(xs,yb)>R*R]+=1
    g=yt>=yb; xs,yt,yb=xs[g],yt[g],yb[g]
    return np.vstack([np.c_[xs,yt],np.c_[xs,yb]])

J2=2*np.pi
print("1) Disque : rotations seules (centre en 0) vs translations aléatoires")
rng=np.random.default_rng(1)
for R in [997,1000,1003,3000,10000,30000]:
    c_rot=(2*np.pi*R-hull_perim(disk_pts(R))[0])*R**(1/3)/J2
    UV=rng.uniform(0,1,(16,2))
    cs=[(2*np.pi*R-hull_perim(disk_pts(R,u,v))[0])*R**(1/3)/J2 for u,v in UV]
    print(f"  R={R:6d}  centre fixe: C={c_rot:.4f}   16 translations: C={np.mean(cs):.4f} ± {np.std(cs)/4:.4f} (err. moy.)")

print("\n2) I_4/3 pour p=1.5 : paramétrisation par x (singularité intégrable)")
def I43_super(p):
    # quart de courbe y=(1-x^p)^(1/p), x in [0,1]; changement x = 1 - s^m pour lisser la singularité
    def f(x):
        y=(1-x**p)**(1/p)
        yp=-x**(p-1)*(1-x**p)**(1/p-1)
        ypp=-(p-1)*x**(p-2)*(1-x**p)**(1/p-2)
        return abs(ypp)**(4/3)/(1+yp*yp)**(2-1/2)  # kappa^{4/3} ds = |y''|^{4/3}/(1+y'^2)^2 * sqrt(1+y'^2) dx
    # intégrer de 0 à 1/2^(1/p) (point symétrique), doubler
    xm=2**(-1/p)
    v,_=quad(f,0,xm,limit=1000,epsabs=1e-13,epsrel=1e-12,points=[1e-6,1e-3])
    return 8*v
for p in [1.5,2.0,2.5,10.0]:
    print(f"  p={p}: I_4/3 = {I43_super(p):.6f}")

print("\n3) R = 30 000, 16 orientations (seed 12345 sur [0, pi/2))")
TH=np.random.default_rng(12345).uniform(0,np.pi/2,16)
R=30000.
for p in [2.5,10.0]:
    Pc=P_cont(R,p); ds=[Pc-hull_perim(hull_points(R,p,t))[0] for t in TH]
    print(f"  superellipse p={p}: C_eff = {np.mean(ds)*R**(1/3)/I43_super(p):.4f} ± {np.std(ds)/4*R**(1/3)/I43_super(p):.4f}")
P,Om,J=invariants(3.0)
ds=[P*R-hull_perim(pts_ellipse(R,3.0,t))[0] for t in TH]
print(f"  ellipse 3:1    : C_eff = {np.mean(ds)*R**(1/3)/J:.4f} ± {np.std(ds)/4*R**(1/3)/J:.4f}")
