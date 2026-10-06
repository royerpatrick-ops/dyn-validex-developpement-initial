"""Test indépendant de la loi ΔP ≈ C·(∫κ^{4/3} ds)·R^{-1/3} et m ≈ c·Ω·R^{2/3}
sur des ellipses (demi-axes R et R/k), moyenne sur 16 orientations."""
import numpy as np  # pyright: ignore[reportMissingImports]
from scipy.spatial import ConvexHull  # pyright: ignore[reportMissingImports]
from scipy.integrate import quad  # pyright: ignore[reportMissingImports]
from balayage_corrige import hull_perim
rng=np.random.default_rng(7); TH=rng.uniform(0,np.pi,16)

def pts_ellipse(R,k,th):
    a,b=R,R/k; c,s=np.cos(th),np.sin(th)
    # ellipse tournée : A x^2 + 2B xy + C y^2 <= 1
    A=c*c/a**2+s*s/b**2; B=c*s*(1/a**2-1/b**2); C=s*s/a**2+c*c/b**2
    xm=np.sqrt(C/(A*C-B*B))
    xs=np.arange(-int(np.ceil(xm))-1,int(np.ceil(xm))+2,dtype=float)
    disc=B*B*xs*xs-C*(A*xs*xs-1); ok=disc>=0; xs,disc=xs[ok],disc[ok]
    y1=(-B*xs+np.sqrt(disc))/C; y0=(-B*xs-np.sqrt(disc))/C
    yt=np.floor(y1+1e-9); yb=np.ceil(y0-1e-9)
    f=lambda x,y:A*x*x+2*B*x*y+C*y*y
    yt[f(xs,yt)>1]-=1; yb[f(xs,yb)>1]+=1
    g=yt>=yb; xs,yt,yb=xs[g],yt[g],yb[g]
    return np.vstack([np.c_[xs,yt],np.c_[xs,yb]])

def invariants(k):  # ellipse unité demi-axes 1, 1/k
    a,b=1.,1./k
    sp=lambda t: np.hypot(a*np.sin(t),b*np.cos(t))
    kap=lambda t: a*b/sp(t)**3
    P=quad(sp,0,2*np.pi,limit=400)[0]
    Om=quad(lambda t:kap(t)**(1/3)*sp(t),0,2*np.pi,limit=400)[0]
    J=quad(lambda t:kap(t)**(4/3)*sp(t),0,2*np.pi,limit=400)[0]
    return P,Om,J

def main():
    print(f"{'a/b':>4} {'R':>6} {'m/(Ω R^2/3)':>12} {'ΔP R^1/3/J':>11} {'D·R^4/3·P/J':>12}")
    for k in [1.5,2.0,3.0,5.0]:
        P,Om,J=invariants(k)
        for R in [300,1000,3000,10000]:
            ms=[];As=[]
            for th in TH:
                Ph,m=hull_perim(pts_ellipse(R,k,th)); ms.append(m); As.append(P*R-Ph)
            mm,AA=np.mean(ms),np.mean(As)
            print(f"{k:4.1f} {R:6d} {mm/(Om*R**(2/3)):12.4f} {AA*R**(1/3)/J:11.4f} {AA/(P*R)*R**(4/3)*P/J:12.4f}")


if __name__ == "__main__":
    main()
