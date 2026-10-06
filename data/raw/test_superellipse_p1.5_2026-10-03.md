# Superellipse p = 1,5 aux grands R — 3 octobre 2026

Script `test_superellipse_singuliere.py`, module `enveloppe_commune.py` corrigé, moyenne sur rotation + translation. Le tableau compare C_eff à C_théo = C∞ − B·(I₃/₂/I₄/₃)·R^(−1/6), avec I₄/₃ = 7,814898 et I₃/₂ = 9,584792.

| R | Tirages | C_eff | C_théo | ΔC | Résidu |
|---|---|---|---|---|---|
| 10³ | 1024 | 0,5937 ± 0,0015 | 0,6070 | −0,0133 | −8,8 σ |
| 3·10³ | 1024 | 0,6170 ± 0,0015 | 0,6257 | −0,0087 | −6,0 σ |
| 10⁴ | 1024 | 0,6348 ± 0,0012 | 0,6427 | −0,0080 | −6,4 σ |
| 3·10⁴ | 1024 | 0,6493 ± 0,0011 | 0,6555 | −0,0063 | −5,6 σ |
| 10⁵ | 1024 | 0,6630 ± 0,0011 | 0,6671 | −0,0041 | −3,9 σ |
| 3·10⁵ | 512 | 0,6734 ± 0,0014 | 0,6758 | −0,0024 | −1,7 σ |
| 10⁶ | 256 | 0,6835 ± 0,0020 | 0,6837 | −0,0002 | −0,1 σ |

## Ajustements

| Modèle pour ΔC | Paramètres | χ² |
|---|---|---|
| Écart nul (loi à deux termes seule) | — | 204,4 (7 points) |
| A·R^(−4/15), exposant prédit | A = −0,084 ± 0,006 | 2,8 (6 ddl) |
| A·R^(−γ), exposant libre | A = −0,076 ± 0,028, γ = 0,254 ± 0,043 | 2,7 (5 ddl) |

## Interprétation

Près des quatre points de courbure infinie, où κ ~ |x|^(p−2), le modèle local cesse de s'appliquer. Il échoue sur une couche de largeur R^((p−1)/(p+1)) mailles, soit R^(−2/(p+1)) en coordonnées réduites. Cette couche apporte au déficit une correction d'ordre R^(−3(p−1)/(p+1)), c'est-à-dire, sur la constante effective :

ΔC ≈ A · R^(−γ), avec γ(p) = 2(4p − 5) / (3(p + 1)).

- **Exposant.** Pour p = 1,5, la prédiction est γ = 4/15 ≈ 0,267 ; la mesure donne 0,254 ± 0,043.
- **Amplitude.** Le coefficient A n'est pas prédit par l'argument.
- **Lien avec la note originale.** Ce terme explique l'écart de −3,3 σ de p = 1,5 dans le tableau des 11 formes : la prédiction est −0,006 vers R ≈ 1,7·10⁴, pour −0,007 ± 0,002 mesuré.

## Conséquence pour l'énoncé

En termes de ΔP, le terme singulier est en R^(−1/3−γ). Il reste négligeable devant le terme en R^(−1/2) si et seulement si γ > 1/6, c'est-à-dire p > 7/5. Plus généralement, pour une singularité κ ~ |s|^(−β), la condition est β < 3/5.

- Pour p = 1,5, la loi à deux termes avec un reste o(R^(−1/2)) reste vraie, mais la convergence est lente.
- Pour 4/3 < p < 7/5, le terme singulier domine le terme en I₃/₂.
- Pour les points de courbure nulle (p > 2), γ > 2/3 : le terme est négligeable.

## Test discriminant à faire

p = 1,75, où la prédiction est γ ≈ 0,485, soit une décroissance bien plus rapide.

## Test p = 1,75 (run de Patrick, 32 processus)

Pour p = 1,75 : I₄/₃ = 6,497336, I₃/₂ = 6,649578, et l'exposant prédit est γ = 0,485.

| R | Tirages | C_eff | C_théo | ΔC | Résidu |
|---|---|---|---|---|---|
| 10³ | 1024 | 0,6238 ± 0,0011 | 0,6255 | −0,0018 | −1,6 σ |
| 3·10³ | 1024 | 0,6410 ± 0,0010 | 0,6412 | −0,0002 | −0,2 σ |
| 10⁴ | 1024 | 0,6534 ± 0,0008 | 0,6554 | −0,0020 | −2,5 σ |
| 3·10⁴ | 1024 | 0,6650 ± 0,0007 | 0,6661 | −0,0011 | −1,7 σ |
| 10⁵ | 1024 | 0,6757 ± 0,0005 | 0,6757 | −0,0000 | −0,1 σ |
| 3·10⁵ | 512 | 0,6831 ± 0,0007 | 0,6830 | +0,0001 | +0,1 σ |
| 10⁶ | 256 | 0,6891 ± 0,0008 | 0,6896 | −0,0005 | −0,6 σ |

Ajustements de ΔC :

| Modèle | Paramètres | χ² |
|---|---|---|
| Écart nul | — | 12,1 (7 points, p ≈ 0,1) |
| Exposant prédit γ = 0,485 | A = −0,060 ± 0,024 | 5,6 (6 ddl) |
| γ = 4/15 | — | 4,7 |
| γ libre | intervalle à 1 σ : [0,11 ; 0,50] | — |

**Lecture.**
- La loi à deux termes tient à ±0,002 près de 10³ à 10⁶ pour cette forme à points de courbure infinie.
- Le terme singulier est à peine détectable, comme l'annonce sa décroissance rapide.
- Ce test est compatible avec la prédiction, mais il ne mesure pas l'exposant.

**Précision théorique.** La contribution de la couche a elle-même une correction relative en R^(−(p−1)/(2(p+1))), soit R^(−0,1) pour p = 1,5. Un ajustement par une seule puissance mesure donc un exposant effectif. Cette correction rend aussi p = 1,4 peu discriminant.

**Test discriminant proposé.** Localiser le déficit par direction normale, avec la formule de Cauchy (fonction d'appui). L'excès doit se concentrer dans une fenêtre angulaire d'environ R^(−2(p−1)/(p+1)) autour des quatre axes, et cette fenêtre doit rétrécir avec R.

## Test de localisation par direction normale (script `localisation_deficit.py`)

**Méthode.** La formule de Cauchy donne ΔP = ∫ [h_RK(φ) − h_P(φ)] dφ. Le déficit mesuré est comparé à la densité locale prédite, C∞·κ_a^(1/3) − B·κ_a^(1/2), angle normal par angle normal. Les intervalles sont logarithmiques en ψ (distance angulaire à l'axe), de 10⁻⁵ à π/4, et les 8 demi-fenêtres sont repliées. Pour chaque tirage, la somme des intervalles redonne P_cont − P_hull à 10⁻⁹ près. Les tirages sont indépendants de ceux du balayage.

**Superellipse p = 1,5, 1024 tirages par rayon.** Excès = déficit mesuré − déficit prédit, cumulé depuis l'axe ; ψ0 = R^(−0,4).

| R | Excès total | Excès à ψ ≤ ψ0/4 | Excès à ψ > 0,1 rad | Mesuré / prédit dans ψ ≤ ψ0/4 |
|---|---|---|---|---|
| 10³ | −0,00752 ± 0,00117 | −0,00971 ± 0,00033 | +0,00129 ± 0,00090 | 0,694 |
| 3·10³ | −0,00483 ± 0,00078 | −0,00550 ± 0,00019 | −0,00005 ± 0,00053 | 0,669 |
| 10⁴ | −0,00376 ± 0,00042 | −0,00304 ± 0,00010 | −0,00063 ± 0,00030 | 0,638 |
| 3·10⁴ | −0,00148 ± 0,00028 | −0,00154 ± 0,00006 | −0,00004 ± 0,00019 | 0,648 |
| 10⁵ | −0,00074 ± 0,00018 | −0,00080 ± 0,00003 | +0,00005 ± 0,00010 | 0,637 |

**Contrôle : ellipse 2:1, 512 tirages.**
- R = 10⁴ : excès total +0,00021 ± 0,00046, et −0,00004 ± 0,00006 près des axes.
- R = 10⁵ : excès total +0,00015 ± 0,00014, et +0,00001 ± 0,00002 près des axes.

**Résultats.**
1. *Localisation.* Tout l'excès se trouve à ψ ≤ ψ0/4 autour des quatre axes. Au-delà de 0,1 rad, l'excès est nul pour R ≥ 3·10³. Le contrôle est nul partout.
2. *Largeur.* La largeur à mi-excès suit ψ_1/2 ∝ R^(−ν) avec ν = 0,378 ± 0,086, pour 0,400 prédit (χ² = 0,7 pour 3 ddl), et ψ_1/2·R^0,4 ≈ 0,04 à 0,055. Les courbes cumulées se superposent en ψ·R^0,4.
3. *Couche autosimilaire.* Dans la fenêtre ψ ≤ ψ0/4, le déficit réel se stabilise à environ 64 % de la prédiction du modèle local (0,694 → 0,637).
4. *Amplitude.* Dans la même fenêtre, l'excès décroît avec un exposant effectif μ = 0,543 ± 0,010 entre 10³ et 10⁵.
   - Une puissance pure R^(−0,6) est rejetée : χ² = 37 pour 4 ddl.
   - La forme R^(−0,6)·(α + β·R^(−0,1)) s'ajuste bien (χ² = 3,4 pour 3 ddl), avec α = 1,13 et β = −1,02. Le terme en R^(−0,1) est celui qu'attend la structure de la couche elle-même, qui comprend la part en R^(−0,7) du terme B.
   - L'exposant asymptotique 0,6 est donc compatible, mais il n'est pas mesuré directement à ces rayons. Le balayage en ΔC donnait μ = 1/3 + γ = 0,587 ± 0,043.

**Conclusion.** Le terme singulier est une couche limite localisée aux points de courbure infinie. Sa largeur angulaire décroît comme R^(−2(p−1)/(p+1)), et le déficit qu'elle porte vaut une fraction fixe (environ 64 %) de ce que prédit le modèle local. Le mécanisme est confirmé ; l'amplitude reste mesurée, non calculée.
