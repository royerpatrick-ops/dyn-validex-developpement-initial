"""Δ(R) = T(R) − C∞, calculé avec la même table de G pour la somme discrète sur Z² et l'intégrale continue."""
import math, numpy as np
from formule_integrale import echantillon
rng=np.random.default_rng(29); g0=4*math.sqrt(2)/15
la=np.linspace(math.log(2e-3),math.log(8),90)
Gt=np.array([echantillon(math.exp(l),600000,rng)[0].mean() for l in la])
def G(a):
    a=np.asarray(a,float); out=np.empty_like(a)
    lo=a<math.exp(la[0]); hi=a>math.exp(la[-1]); mid=~lo&~hi
    out[lo]=g0*a[lo]**-1.5; out[hi]=(1/9)*a[hi]**-3*(1+1.8/a[hi]**2)
    out[mid]=np.exp(np.interp(np.log(a[mid]),la,np.log(np.maximum(Gt,1e-300))))
    return out
# intégrale continue (6/π²)∫ a G(a) da, mêmes formes
u=np.linspace(math.log(1e-9),math.log(1e4),400001); a=np.exp(u)
Cinf=6/math.pi**2*np.trapezoid(a*a*G(a),u) if hasattr(np,"trapezoid") else 6/math.pi**2*np.trapz(a*a*G(a),u)
print(f"C∞ de la table : {Cinf:.4f}")
def T(R,amax=200):
    eps=R**(-1/3); X=int(amax/eps)+1; tot=0.0; p=np.arange(-X,X+1)
    for q in range(0,X+1):
        pp=p if q>0 else p[p>0]; r2=pp*pp+q*q; m=(r2<=X*X)&(np.gcd(pp,q)==1)
        tot+=2*np.sum(G(eps*np.sqrt(r2[m].astype(float))))
    tail=6/math.pi**2*(1/9)/amax
    return R**(-2/3)/(2*math.pi)*tot+tail
B=g0*3.857623/(2*math.pi)
print("     R      x=R^-1/6   Δ(R)      −B·x      Δ+B·x")
for R in [1000,3000,10000,30000,100000,300000,1000000]:
    d=T(R)-Cinf; x=R**(-1/6); print(f"  {R:>8d}  {x:.4f}  {d:+.5f}  {-B*x:+.5f}  {d+B*x:+.5f}",flush=True)
