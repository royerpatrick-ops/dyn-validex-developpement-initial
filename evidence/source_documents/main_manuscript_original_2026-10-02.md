# Déficit de périmètre de l'enveloppe entière d'un convexe lisse : une loi en ∫κ^{4/3}

Oct 2, 2026 · @Patrick

## Résumé

Soit K un convexe plan à bord lisse et P\_R l'enveloppe convexe des points de ℤ² contenus dans RK. Nous conjecturons que le déficit de périmètre vérifie ΔP(R) ≈ C · (∫κ^{4/3} ds) · R^{−1/3}, où l'intégrale porte sur le bord de K et où C ne dépend pas de la forme. Sa valeur limite est C∞ = 0,719 ± 0,001, donnée par une intégrale explicite, et elle est approchée lentement, avec une correction universelle en R^{−1/6} :

```latex
\Delta P(R) \simeq C_\infty \Bigl(\int_{\partial K}\kappa^{4/3}\,ds\Bigr) R^{-1/3} \;-\; B \Bigl(\int_{\partial K}\kappa^{3/2}\,ds\Bigr) R^{-1/2}, \qquad C_\infty = 0{,}7192,\; B = 0{,}2894
```

Les deux constantes sont calculées, et non ajustées. Le premier terme se déduit d'un modèle local : une transformation affine qui conserve les aires ramène le bord à une parabole et ℤ² à un réseau aléatoire de mesure de Haar. Réécrit comme une somme sur les vecteurs primitifs du réseau, ce modèle donne C∞ sous forme d'intégrale explicite, dont les deux comportements limites sont exacts (4√2/15·a^{−3/2} et 1/(9a³)) : **C∞ = 0,7192 ± 0,0003**. Il redonne au passage la constante de Balog et Deshouillers pour le nombre de sommets (0,54971 contre 0,54967). Le second terme vient de l'arithmétique de ℤ², dont les vecteurs entiers ne sont jamais plus courts que 1. Il a deux sources, toutes deux calculables : la fonction zêta d'Epstein de ℤ² pour les arêtes courtes (0,2315), et la discrétude des vecteurs partenaires pour les arêtes à deux colonnes (0,0579), soit **B = 0,2894**.

La loi est confrontée à des calculs exacts d'enveloppes : 11 formes pour R de 100 à 30 000 (64 tirages par point), puis 4 formes jusqu'à R = 10⁶ (128 tirages par point). Sans aucun paramètre ajusté, elle reproduit les 92 mesures de R = 1 230 à 10⁶ avec χ²/n = 1,16. Tous les calculs ont été reproduits à l'identique sur une seconde machine. Nous n'avons trouvé cette loi ni dans les travaux sur l'enveloppe entière, ni en géométrie discrète.

## 1. Introduction et contexte

L'enveloppe entière P\_R = conv(RK ∩ ℤ²) est étudiée depuis trente ans, mais presque uniquement à travers le nombre de ses sommets et le déficit d'aire. Sur son périmètre, la littérature ne donne à notre connaissance ni ordre de grandeur démontré, ni dépendance explicite en la forme.

Ce qui est établi :

- **Ordre R^{2/3} du nombre de sommets.** Balog et Bárány l'ont montré pour le disque ([Cowles DP 962, 1990](https://elischolar.library.yale.edu/cowles-discussion-paper-series/1205/)). Bárány et Larman l'ont étendu aux corps à courbure bornée en dimension d, avec l'ordre r^{d(d−1)/(d+1)} ([Math. Ann., 1998](https://link.springer.com/article/10.1007/s002080050217) ; [synthèse](https://users.renyi.hu/~barany/cikkek/psurveyAMS.pdf)).
- **Constante moyenne pour le disque.** Balog et Deshouillers ont obtenu une moyenne sur r de m·r^{−2/3} tendant vers 3,453… ([rappel](https://arxiv.org/html/2602.06897)).
- **Courbe lisse générale.** Huxley a montré que la moyenne sur les translations vaut B·R^{2/3}, avec B explicite ([Period. Math. Hungar. 68, 2014](https://link.springer.com/article/10.1007/s10998-014-0024-5)). Nous n'avons pas pu consulter la forme de B.
- **Géométrie discrète.** Pour une forme lisse numérisée au pas h, les sommets de l'enveloppe sont en Θ(h^{−2/3}) et à distance O(h^{3/2}) du bord ([Lachaud, Coeurjolly, Roussillon, DGMM 2025](https://perso.liris.cnrs.fr/david.coeurjolly/publication/dgmm-25-geometryconvex/dgmm-25-geometryconvex.pdf)). Les estimateurs de longueur MLP et λ-MST convergent en O(h^{4/3}) dans les expériences, pour des bornes démontrées en O(h) ou moins bonnes ([Roussillon et al., 2012](https://perso.liris.cnrs.fr/tristan.roussillon/static/file/2012bookchapter.pdf)).

Cette note propose une loi asymptotique explicite pour le périmètre et la teste numériquement. Elle vérifie aussi, au passage, la proportionnalité du nombre de sommets au périmètre affine.

## 2. Énoncé

Soit K un convexe compact du plan dont le bord est de classe C² par morceaux, et dont la courbure κ ne s'annule ou ne diverge qu'en un nombre fini de points, avec ∫κ^{4/3} ds fini. Posons P\_R = conv(RK ∩ ℤ²), avec une position générique : une rotation ou une translation aléatoire de K, en moyenne.

**Conjecture.** Quand R → ∞, en moyenne sur la position de K par rapport au réseau (rotations et translations), le déficit de périmètre vérifie :

```latex
P(RK) - P(P_R) \;=\; C_\infty \Bigl(\int_{\partial K} \kappa^{4/3}\, ds\Bigr) R^{-1/3} \;-\; B \Bigl(\int_{\partial K} \kappa^{3/2}\, ds\Bigr) R^{-1/2} \;+\; o\bigl(R^{-1/2}\bigr)
```

Ici C∞ et B sont des constantes universelles, indépendantes de K, toutes deux définies par des formules explicites : C∞ = 0,7192 ± 0,0003 par l'intégrale de la §3.4, et B = 0,2894 par le calcul arithmétique de la §3.5. Les simulations, ajustées librement, donnaient C∞ = 0,7204 ± 0,0005 et B = 0,296 ± 0,003 ; avec C∞ fixé, B = 0,289 ± 0,001. Le déficit relatif est donc en R^{−4/3}, et la constante effective C\_eff(R) = ΔP·R^{1/3} / ∫κ^{4/3} ds tend vers C∞ par valeurs inférieures, avec un écart en R^{−1/6}.

**Loi compagne (nombre de sommets).** Dans les mêmes conditions :

```latex
f_0(P_R) \;\sim\; \frac{c_0}{2\pi}\,\Omega(K)\, R^{2/3}, \qquad \Omega(K) = \int_{\partial K} \kappa^{1/3}\, ds
```

Ici c₀ ≈ 3,453 est la constante de Balog et Deshouillers pour le disque, et Ω(K) est le périmètre affine. Cette loi est probablement contenue dans le résultat de Huxley (2014), à confirmer sur le texte. Le modèle local de la §3.2 la reproduit sans aucun paramètre ajusté.

**Cas particuliers.** Le disque donne ∫κ^{4/3} ds = ∫κ^{3/2} ds = 2π, donc ΔP ≈ 4,52 · R^{−1/3} − 1,82 · R^{−1/2}. Les superellipses avec p > 2 ont des points de courbure nulle, et celles avec 1 < p < 2 des points de courbure infinie ; les deux obéissent à la même loi, sans changement d'exposant, tant que les intégrales de courbure sont finies. Un polygone relève d'un autre régime (ΔP = O(1), soit un déficit relatif en R^{−1}), car ∫κ^{4/3} y diverge au voisinage des coins.

## 3. Origine de la loi

### 3.1 Argument heuristique : la calotte d'aire unité

Toutes les échelles locales découlent d'une seule condition : une calotte du bord contient de l'ordre d'un point du réseau quand son aire est de l'ordre de la maille. Soit κ\_a = κ/R la courbure réelle au point considéré. Une calotte de largeur w a une profondeur δ ≈ κ\_a w²/8 et une aire de l'ordre de κ\_a w³. L'égaler à 1 fixe l'échelle locale :

```latex
w \sim \kappa_a^{-1/3}, \qquad \delta \sim \kappa_a^{1/3}
```

De là découlent les deux lois :

1. **Sommets.** La densité de sommets par unité de longueur est ρ ∝ 1/w ∝ κ\_a^{1/3}. En intégrant sur le bord de RK, on obtient f₀ ∝ ∫κ\_a^{1/3} ds\_a = Ω(K) · R^{2/3}.
2. **Périmètre.** Deux effets, du même ordre, retranchent du périmètre :
   - *corde contre arc* : chaque arête de longueur w perd κ\_a² w³ / 24, soit κ\_a² w² / 24 ∝ κ\_a^{4/3} par unité de longueur ;
   - *retrait vers l'intérieur* : les sommets sont à une profondeur δ ∝ κ\_a^{1/3} du bord, et un décalage δ vers l'intérieur réduit la longueur de κ\_a δ ∝ κ\_a^{4/3} par unité de longueur.

En intégrant, ΔP ∝ ∫κ\_a^{4/3} ds\_a = R^{−1/3} ∫κ^{4/3} ds. L'exposant et l'invariant de forme sont donc fixés par cette seule échelle locale, quel que soit le poids relatif des deux effets.

Cet argument fixe les échelles, mais pas la constante : un espacement régulier des sommets donnerait C ≈ 0,14, soit cinq fois moins que la valeur mesurée. La constante dépend de la loi jointe des longueurs d'arêtes et des retraits, que fournit le modèle suivant.

### 3.2 Modèle local : un réseau aléatoire sous une parabole

Au voisinage d'un point du bord de courbure réelle κ\_a = κ/R, l'application x' = κ\_a^{1/3} x\_t, y' = κ\_a^{−1/3} y\_n (coordonnées tangente et normale) conserve les aires. Elle envoie le bord sur la parabole y = −x²/2 et ℤ² sur un réseau d'aire 1. Le cisaillement qui fait glisser la parabole le long d'elle-même agit sur ce réseau comme le flot horocyclique, et l'augmentation de R comme le flot géodésique. Par équirépartition des longs horocycles, le réseau local se comporte en moyenne comme un réseau aléatoire distribué selon la mesure de Haar, translaté aléatoirement.

Dans ces coordonnées, une courbe intérieure à distance δ(s) du bord est plus courte de ∫(κδ − δ'²/2) ds au second ordre. Le déficit par unité de longueur réelle vaut donc κ\_a^{4/3} · D, où D = ⟨g − g'²/2⟩ et g est l'écart entre la parabole et l'enveloppe du réseau aléatoire. En intégrant le long du bord, on obtient exactement le premier terme de la conjecture, avec **C∞ = D**, une espérance de Haar. De même, la densité de sommets par unité de longueur normalisée doit valoir c₀/2π.

Nous avons calculé ces deux espérances par Monte-Carlo (§5.4). Le modèle redonne c₀/2π à 0,02 % près et fournit C∞ = 0,718 ± 0,003, sans paramètre ajusté sur les simulations du problème réel.

### 3.3 D'où vient le terme correctif

Le modèle local suppose un réseau aléatoire. Le vrai réseau ℤ² en diffère sur un point : ses vecteurs sont entiers, donc de longueur au moins 1, et leurs longueurs ne prennent que des valeurs discrètes (1, √2, 2, √5…). Vu à travers la transformation locale, qui contracte la direction tangente d'un facteur ε = κ\_a^{1/3} = (κ/R)^{1/3}, cela affecte les arêtes à deux extrémités :

- **les arêtes de petit a**, portées par un vecteur w₀ ∈ ℤ² court et presque tangent : de nombreux points alignés sur une même droite. Leur nombre est fixé par la somme discrète sur les vecteurs primitifs de ℤ², et non par l'intégrale continue ;
- **les arêtes de grand a**, à deux colonnes : elles exigent un vecteur partenaire v (det(w₀, v) = 1) presque perpendiculaire à w₀, ce qui est impossible si |v| < 1. Elles disparaissent au-delà de a ≈ √2·(R/κ)^{1/6}.

Les deux effets sont du même ordre, ε^{1/2} = (κ/R)^{1/6}. Intégrés le long du bord avec le poids κ\_a^{4/3}, ils donnent le second terme de la conjecture, en R^{−1/2} ∫κ^{3/2} ds. Leurs coefficients sont calculés en §3.5.

### 3.4 C∞ comme intégrale explicite

Le modèle local se réécrit comme une somme sur les vecteurs primitifs du réseau, à la manière de Huxley pour le nombre de sommets. Pour chaque vecteur primitif w, les points du réseau se rangent sur des droites parallèles à w. Le côté de l'enveloppe de direction w se trouve sur la première de ces droites, en partant du haut, qui contient un point sous la parabole : c'est une arête si elle en contient au moins deux, un simple sommet sinon. Cette description est exacte.

Un cisaillement qui préserve la parabole rend w horizontal, de longueur a. Les droites sont alors horizontales, espacées de 1/a, et la configuration est décrite par trois paramètres uniformes :

- δ ∈ \[0, 1/a) : la profondeur de la première droite sous le sommet de la parabole ;
- s ∈ \[0, a) : la position des points sur cette droite ;
- β ∈ \[0, a) : le décalage des points d'une droite à la suivante.

La droite numéro k est à la profondeur d\_k = δ + k/a et porte les points x = s + kβ + ja, j ∈ ℤ. Si la première droite qui coupe la parabole y porte des points de x₁ à x₂ (x₁ < x₂), l'arête correspondante apporte au déficit :

```latex
c(a,\delta,s,\beta) \;=\; d\,(x_2 - x_1) \;-\; \frac{x_2^3 - x_1^3}{3}
```

Le théorème de la moyenne de Siegel, appliqué aux couples (réseau, vecteur primitif), transforme la somme sur tous les w en intégrale. Le facteur 6/π² = 1/ζ(2) est la probabilité que deux entiers soient premiers entre eux :

```latex
C_\infty \;=\; \frac{6}{\pi^2}\int_0^{\infty} a\; \mathbb{E}_{\delta,s,\beta}\bigl[c(a,\delta,s,\beta)\bigr]\, da, \qquad \frac{c_0}{2\pi} \;=\; \frac{6}{\pi^2}\int_0^{\infty} a\; \mathbb{P}_{\delta,s,\beta}\bigl[\text{arête}\bigr]\, da
```

La seconde formule, pour le nombre de sommets, sert de contrôle. L'intégrande a deux comportements limites connus :

- **a petit** (points serrés sur la droite) : la première droite est toujours atteinte, et l'on obtient exactement E\[c\] ≈ g₀·a^{−3/2} avec **g₀ = 4√2/15** ;
- **a grand** : seules contribuent les configurations à deux colonnes de points à ±a/2. En réduisant à l'échelle 1/a autour de la profondeur a²/8, la condition d'arête devient, pour chaque inclinaison m ∈ (−1, 1) des colonnes, deux triangles d'aire totale (1 − m²)/2. Leur intégrale vaut J = 2/3 ; chaque arête apporte a³/24 et la probabilité d'une arête vaut 4J/a⁶. D'où exactement **E\[c\] ≈ 1/(9a³)**, avec une correction relative en 1/a².

La formule se simplifie donc en :

```latex
C_\infty \;=\; \frac{6}{\pi^2}\Bigl[\int_0^{A} a\,\mathbb{E}[c](a)\,da \;+\; \frac{1}{9A}\Bigr] \;+\; O(A^{-3})
```

où seule la zone intermédiaire, a ≈ 0,3 à 3, demande un calcul numérique ; on ne lui connaît pas de forme close, comme pour la constante c₀ de Balog et Deshouillers. L'évaluation (§5.4) donne c₀/2π = 0,54971 ± 0,00012, contre 0,54967 attendu, et **C∞ = 0,7192 ± 0,0003**. À l'intérieur du modèle local, la dérivation est rigoureuse : la description géométrique est exacte et le théorème de Siegel est classique.

### 3.5 Calcul du terme correctif B

Pour ℤ², une arête de direction primitive w₀ n'apparaît qu'au point du bord où la tangente est parallèle à w₀, avec a = ε|w₀|. La constante effective s'écrit donc exactement comme une somme discrète :

```latex
C_{\mathrm{eff}}(R) \;=\; \frac{\varepsilon^2}{2\pi}\sum_{w_0\in\mathbb{Z}^2\ \text{primitif}} \mathbb{E}[c]\bigl(\varepsilon |w_0|\bigr), \qquad \varepsilon = \kappa_a^{1/3}
```

que l'intégrale continue (6/π²)∫d²w₀ transforme en C∞. Le terme correctif est l'écart entre cette somme et cette intégrale. Il a deux sources.

**Petits a.** Avec E\[c\] ≈ g₀ a^{−3/2}, l'écart entre la somme sur les vecteurs primitifs de |w₀|^{−3/2} et son intégrale est la valeur régularisée de la fonction zêta d'Epstein de ℤ² restreinte aux vecteurs primitifs :

```latex
B_{\text{petits}} \;=\; \frac{g_0}{2\pi}\,\bigl|Z_{\mathrm{prim}}(3/2)\bigr|, \qquad Z_{\mathrm{prim}}(3/2) = \frac{4\,\zeta(3/4)\,\beta(3/4)}{\zeta(3/2)} = -3{,}8576
```

où β est la fonction bêta de Dirichlet. On obtient **B\_petits = 0,2315**.

**Grands a.** Une arête à deux colonnes exige |m| < 1 avec m = λ·|v|, où λ = a²ε/2 et v est le vecteur partenaire entier (det(w₀, v) = 1). Comme |v| ne prend que des valeurs discrètes ≥ 1, la probabilité d'arête est multipliée par le rapport

```latex
\rho(\lambda) \;=\; \frac{\pi\lambda}{8}\sum_{\substack{v\ \text{primitif}\\ |v|<1/\lambda}} \frac{1-\lambda^2|v|^2}{|v|}, \qquad B_{\text{grands}} \;=\; -\frac{6}{\pi^2}\cdot\frac{1}{9}\cdot\frac{\sqrt{2}}{4}\int_0^{\infty}\lambda^{-3/2}\bigl(\rho(\lambda)-1\bigr)\,d\lambda
```

qui tend vers 1 pour λ petit et s'annule pour λ ≥ 1. L'intégrale vaut −2,4238, d'où **B\_grands = 0,0579**.

**Total : B = 0,2315 + 0,0579 = 0,2894**, à comparer à 0,289 ± 0,001 ajusté sur les simulations. La somme discrète a aussi été vérifiée classe par classe sur le disque réel : pour chaque longueur |w₀| de 1 à √34, la contribution mesurée suit la prédiction du régime des petits a à quelques pour cent près, à R = 10⁵ comme à 10⁶.

### 3.6 Une identité exacte : le théorème de Pick

Pour un polygone à sommets entiers, le théorème de Pick donne aire = I + B/2 − 1, où I et B comptent les points du réseau intérieurs et sur le bord. Moyenné sur les translations, le nombre total de points du réseau contenus dans RK vaut exactement l'aire de RK. On en déduit, **exactement et pour tout R** :

```latex
\mathbb{E}\bigl[\mathrm{aire}(RK) - \mathrm{aire}(P_R)\bigr] \;=\; \frac{\mathbb{E}[B]}{2} + 1
```

Le déficit d'aire se ramène donc au nombre de points du réseau situés sur le bord de l'enveloppe. Comme la densité de déficit de périmètre s'écrit D = ⟨g⟩ − ½⟨g'²⟩, où ⟨g⟩ est précisément le déficit d'aire par unité de longueur, la constante se décompose en :

```latex
C_\infty \;=\; \frac{n_\partial}{2} \;-\; \frac{1}{2}\,\langle g'^2\rangle
```

où n\_∂ est la densité, par unité de longueur normalisée, des points du réseau sur le bord de l'enveloppe. Dans le modèle local, on mesure ⟨g⟩ = 1,1509 ± 0,0008 directement et n\_∂/2 = 1,1499 ± 0,0002 par comptage, en accord avec l'identité, d'où ½⟨g'²⟩ ≈ 0,432. Le bord porte en moyenne environ 4,2 points du réseau par sommet, à cause des arêtes qui alignent plusieurs points.

**Statut de C∞.** La partie intermédiaire de l'intégrale (a ≈ 0,3 à 3) repose sur la probabilité qu'une calotte ne contienne aucun point du réseau. Pour une direction fixée, la mesure de l'ensemble des paramètres β qui conviennent est polynomiale par morceaux, avec des raccords aux fractions de Farey. C'est la structure qui définit la constante c₀ de Balog et Deshouillers, connue seulement comme une intégrale sur un domaine polygonal. C∞ a vraisemblablement le même statut : une intégrale explicite, sans forme close élémentaire connue.

## 4. Méthode numérique

Chaque enveloppe est calculée exactement, sans grille bornée. Pour chaque abscisse entière x, on détermine l'intervalle des ordonnées entières contenues dans RK, puis on prend l'enveloppe convexe des deux extrémités de chaque colonne. Le coût est en O(R log R).

- **Colonnes.** Pour chaque abscisse du réseau décalé ℤ² − (u, v), on calcule les bornes de RK. Pour une superellipse tournée, le minimum de F(x, ·) est trouvé par section dorée et les bornes par bisection (75 itérations). Pour une ellipse, elles sont données en forme fermée. Les bornes sont ensuite corrigées à l'entier par un test exact. Sur 33 cas test (11 formes, 3 rayons, rotation et translation aléatoires), le résultat est identique à l'enveloppe calculée sur une grille pleine.
- **Périmètre continu et invariants.** P, Ω = ∫κ^{1/3} ds et I₄/₃ = ∫κ^{4/3} ds sont calculés par quadrature adaptative. Pour les superellipses, on paramètre par x avec le changement x = s⁶, qui supprime les singularités de courbure. On obtient ainsi I₄/₃ = 7,8149 pour p = 1,5 (et non 7,768 comme dans la première version). Contrôle : I₄/₃ = Ω = P = 2π pour le disque, à 10⁻¹² près.
- **Moyennes.** Pour chaque couple (forme, R), on tire 64 échantillons. Chacun combine une rotation θ ∈ \[0, π/2) et un décalage (u, v) ∈ \[0, 1)² du réseau. Les graines sont déterministes, dérivées de (forme, paramètre, R). Chaque C\_eff est donné avec son erreur sur la moyenne.
- **Jeu de formes.** Il comprend 7 superellipses, avec p ∈ {1,5 ; 2 ; 2,5 ; 3 ; 4 ; 5 ; 10}, et 4 ellipses, de rapports d'axes 1,5 ; 2 ; 3 et 5. On utilise 26 rayons, répartis géométriquement de 100 à 30 000.

**Mise en garde méthodologique.** Une version antérieure de ces calculs rognait la forme sur le carré \[−R−1, R+1\]². Pour p > 2, cela crée un faux plateau du déficit relatif. Tout générateur de masque doit vérifier que la forme ne touche pas le bord de la grille. Un second piège concerne la première version de cette note : la moyenne y portait sur les seules rotations autour d'un centre entier. Pour le disque, une rotation laisse l'ensemble des points du réseau inchangé, si bien que les 16 tirages étaient identiques et que la valeur publiée n'était qu'un tirage unique. Il en résultait par exemple un faux creux à R = 1 000 (C\_eff = 0,468, contre 0,637 ± 0,005 avec des translations). D'où la moyenne sur les translations adoptée ici.

**Compléments pour les grands R et le modèle local.**

- **Grands R.** Le disque, les ellipses 2:1 et 5:1 et la superellipse p = 3 sont calculés pour 9 rayons de 30 000 à 10⁶, avec 128 tirages (rotation et translation) par point. Les tirages sont découpés en paquets de 16 pour le calcul parallèle, avec reprise après interruption.
- **Décomposition par arêtes.** Pour le disque, chaque arête apporte exactement RΔφ − |e| au déficit, où Δφ est l'angle sous lequel on la voit du centre. On classe ces contributions selon la longueur normalisée ℓ' = |e|·R^{−1/3}.
- **Modèle local.** Le réseau est tiré selon la mesure de Haar sur SL₂(ℝ)/SL₂(ℤ) : un point z = x + iy du domaine fondamental avec la densité dx dy / y², une rotation uniforme et une translation uniforme. Le tirage est stratifié en 25 tranches de y, de √3/2 à 2²⁴. Seuls les points extrêmes de chaque droite du réseau sous la parabole sont retenus, par une résolution stable de l'équation du second degré suivie d'une vérification exacte. D est estimé en attribuant à chaque arête sa contribution totale en son milieu, sur une fenêtre de longueur 1. L'enveloppe est identique à une énumération brute sur 150 réseaux, dont 100 aux orientations critiques, et le résultat ne dépend pas de la largeur de fenêtre.

## 5. Résultats

Les exposants sont universels, et la forme n'agit que par deux intégrales de courbure : ∫κ^{1/3} ds pour les sommets, ∫κ^{4/3} ds pour le périmètre. La constante du périmètre converge lentement vers C∞ ≈ 0,72, avec une correction en R^{−1/6} dont l'amplitude est elle aussi fixée par une intégrale de courbure.

### 5.1 Exposants et nombre de sommets

Pour les 11 formes, le déficit relatif décroît avec une pente γ\_eff comprise entre 1,305 et 1,314 sur \[300, 30 000\], légèrement sous 4/3 en raison de la dérive de C\_eff (§5.3), y compris pour les points de courbure nulle (p > 2) ou infinie (p = 1,5). Aucun plateau n'apparaît jusqu'à R = 30 000.

Pour R ≥ 1 000, le rapport f₀ / (Ω(K) · R^{2/3}) vaut 0,5502 en moyenne, entre 0,549 et 0,553 selon la forme, contre c₀/2π = 0,5496. Pour le disque en translation, on obtient 0,5498, en accord à 0,04 % avec la constante de Balog et Deshouillers.

### 5.2 Quel invariant de forme ?

&#91;embedded content: Calcul de l'auteur · 7 superellipses et 4 ellipses, 64 tirages (rotation + translation) chacune, R de 10 000 à 30 000\]

Parmi les invariants ∫κ^a ds, seul a = 4/3 rend C\_eff presque identique d'une forme à l'autre : 1,2 % de dispersion pour R ≥ 10 000, contre 4,3 % pour a = 5/4 et 6,5 % pour a = 17/12. Le tableau donne la constante forme par forme, avec l'erreur sur la moyenne.

| Forme | ∫κ^{4/3} ds | C\_eff sur \[10 000, 30 000\] |
| --- | --- | --- |
| Disque (p = 2) | 6,2832 | 0,6623 ± 0,0009 |
| Superellipse p = 2,5 | 6,4047 | 0,6618 ± 0,0012 |
| Superellipse p = 3 | 6,6735 | 0,6590 ± 0,0013 |
| Ellipse a/b = 1,5 | 6,9299 | 0,6586 ± 0,0013 |
| Ellipse a/b = 2 | 7,6897 | 0,6572 ± 0,0014 |
| Superellipse p = 4 | 7,2373 | 0,6567 ± 0,0014 |
| Superellipse p = 5 | 7,7517 | 0,6532 ± 0,0014 |
| Superellipse p = 10 | 9,6942 | 0,6494 ± 0,0016 |
| Ellipse a/b = 3 | 9,2667 | 0,6492 ± 0,0016 |
| Superellipse p = 1,5 | 7,8149 | 0,6431 ± 0,0022 |
| Ellipse a/b = 5 | 12,2867 | 0,6366 ± 0,0019 |

À R fini, les écarts entre formes dépassent l'erreur statistique (dispersion de 1,2 % contre un bruit de 0,2 %) : les formes dont la courbure varie le plus (ellipse 5:1, p = 1,5, p = 10) sont en dessous. Ces écarts ne sont pas une violation de l'universalité : ils viennent du terme correctif en ∫κ^{3/2} ds, qui pèse davantage pour ces formes, et ils disparaissent quand on en tient compte (§5.3).

### 5.3 Dérive de C\_eff avec R

&#91;embedded content: Calcul de l'auteur · 11 formes, 14 rayons affichés sur 26 calculés, 64 tirages (rotation + translation) chacun\]

C\_eff n'a pas atteint de plateau à R = 30 000 : la moyenne des 11 formes passe de 0,575 sur \[100, 300\] à 0,653 sur \[10 000, 30 000\], et continue de monter jusqu'à 0,676 – 0,691 à R = 10⁶. Plusieurs lois de correction ont été comparées sur les 4 formes calculées jusqu'à 10⁶ (R ≥ 1 000) :

| Loi de correction | χ²/ddl | Verdict |
| --- | --- | --- |
| C∞ − b·R^{−1/3} | 3,19 | rejetée : sa prédiction à l'aveugle depuis R ≤ 30 000 est trop basse de 4 à 7 σ à R = 10⁶ |
| a + b·ln R (sans limite) | 3,48 | rejetée : la pente en ln R décroît d'un facteur 2 entre 10³ et 10⁶ |
| C∞ + b/ln R | 1,18 | compatible, mais sans mécanisme |
| C∞ − b·R^{−γ}, γ libre | 1,09 | γ = 0,172 ± 0,012, compatible avec 1/6 |
| **C∞ − b·R^{−1/6}** | **1,08** | **prédite par le mécanisme de la §3.3** |

Avec γ = 1/6, l'amplitude b varie selon la forme exactement comme ∫κ^{3/2} ds / ∫κ^{4/3} ds : rapportée à cet invariant, elle vaut 0,288 ; 0,292 ; 0,287 et 0,292 pour le disque, l'ellipse 2:1, l'ellipse 5:1 et la superellipse p = 3. Les 4 formes suivent alors la droite calculée aux §3.4 et 3.5, sans aucun paramètre ajusté :

&#91;embedded content: Calcul de l'auteur · R ≤ 30 000 : 64 tirages par point ; R > 30 000 : 128 tirages par point, répliqués sur une seconde machine\]

**Test sans paramètre libre.** Avec C∞ = 0,7192 et B = 0,2894, tous deux calculés, la loi reproduit les 92 points de mesure (4 formes, R de 1 230 à 10⁶) avec χ² = 106,8, soit χ²/n = 1,16. Aucun point ne s'écarte de plus de 2,7 σ. Pour comparaison, l'ajustement libre des deux constantes donne C∞ = 0,7204 ± 0,0005 et B = 0,296 ± 0,003 (χ²/ddl = 1,13) : il n'apporte presque rien, ce qui confirme que les deux constantes calculées sont les bonnes.

Un contrôle plus fin porte sur la somme discrète de la §3.5. En classant les arêtes du disque réel par la longueur primitive a = ε|w₀|, la contribution mesurée coïncide tranche par tranche avec la prédiction du modèle sommé sur les vrais vecteurs de ℤ² (total 0,6765 des deux côtés à R = 10⁵). Le seul écart notable concerne a > 8,5 : le disque réel n'y a presque plus d'arêtes, comme le prévoit la coupure à a ≈ √2·R^{1/6} ≈ 9,6 due aux vecteurs partenaires entiers.

### 5.4 Modèle local

Le modèle de la §3.2 a été évalué de deux façons indépendantes, chacune reproduite à l'identique sur deux machines : par Monte-Carlo direct (25 tranches de profondeur y, 20 000 réseaux par tranche), et par la formule intégrale de la §3.4. Pour celle-ci, on intègre directement pour a ≤ 6, avec 2 × 10⁶ tirages par nœud et une queue analytique en a → 0, et on échantillonne les configurations à deux colonnes pour a > 6.

| Quantité | Monte-Carlo direct | Formule intégrale | Référence |
| --- | --- | --- | --- |
| Sommets par unité de longueur normalisée | 0,5503 ± 0,0019 | 0,54971 ± 0,00012 | c₀/2π = 0,54967 (Balog et Deshouillers) |
| Queue des grands a | part en y^{−0,2500 ± 0,0015} | a³·E\[c\] = 0,1111 ± 0,0005 à a = 16 | 1/9 = 0,1111 (calcul exact, §3.4) |
| C∞ | 0,718 ± 0,003 | **0,7192 ± 0,0003** | 0,7204 ± 0,0005 (ajustement libre des simulations) |

Dans la formule intégrale, les configurations avec a < e^{−12} apportent 0,0011, celles avec e^{−12} < a < 6 apportent 0,7066, et celles avec a > 6 apportent 0,0115. La valeur a³·E\[c\] vaut 0,1170 à a = 6 et 0,1138 à a = 8, puis se stabilise à 1/9 : la correction est en 1/a², comme prévu. En tenant compte des erreurs systématiques (discrétisation de l'intégrale, raccord à a = 6), on retient C∞ = 0,7192 ± 0,0003.

Le calcul de B (§3.5) a été évalué de la même façon : la régularisation de la somme sur les vecteurs primitifs de ℤ² est contrôlée numériquement (somme jusqu'à |w| = 1 500 : −3,8580, contre −3,8576 attendu), et la somme sur les vecteurs partenaires est faite jusqu'à |v| = 4 000. L'écart Δ(R) entre somme discrète et intégrale, calculé avec la même table de E\[c\], suit −B\_petits·R^{−1/6} à 0,0003 près dès R = 10⁴.

### 5.5 Décomposition par longueur d'arête

Pour le disque, la part de C\_eff portée par les arêtes de longueur normalisée ℓ' < 8 est la même à R = 10⁴, 10⁵ et 10⁶ (environ 0,65). Toute la croissance de C\_eff vient des arêtes plus longues, dont la longueur maximale suit ℓ'max = 2,6 · R^{1/6} à 2 % près de R = 10³ à 10⁶ (2√2 ≈ 2,8 prédit). Au-delà de ℓ' ≈ 4, la part portée par chaque demi-octave de longueur est divisée par un facteur qui tend vers √2, soit une queue en 1/ℓ', l'équivalent de y^{−1/4}. Ces arêtes très longues portent de nombreux points alignés : leur vecteur primitif w₀ est court (a petit), et leur longueur vaut au plus 2√(2/a). La condition a ≥ ε = R^{−1/3} (aucun vecteur entier plus court que 1) donne exactement la borne ℓ'max = 2√2·R^{1/6} : c'est l'effet des petits a de la §3.5, vu sous l'angle des longueurs d'arêtes.

## 6. Discussion, limites et questions ouvertes

La loi est désormais établie dans ses deux termes, et chaque constante est calculée à partir d'un mécanisme identifié. C∞ est l'espérance de Haar d'un modèle local, écrite comme une intégrale explicite ; B est l'écart entre la somme sur les vrais vecteurs de ℤ² et cette intégrale, qui fait intervenir la fonction zêta d'Epstein de ℤ². Ce qui manque, ce sont des démonstrations rigoureuses du passage au modèle local.

**Limites.**

- **Moyenne sur la position.** La loi porte sur la moyenne sur les rotations et translations de K, à l'image du résultat de Huxley pour les sommets. Pour une position fixe, ΔP fluctue fortement (fluctuations diophantiennes).
- **Constantes numériques.** C∞ = 0,719 est défini par une intégrale explicite (§3.4), évaluée numériquement ; on ne lui connaît pas de forme close. B = 0,2894 est calculé (§3.5) ; sa contribution des grands a suppose que la projection du vecteur partenaire est répartie uniformément en dehors des valeurs discrètes imposées par ℤ², ce qui est plausible et confirmé par les données, mais non démontré.
- **Classe de formes.** Superellipses et ellipses seulement, toutes à symétrie centrale. Il manque des formes sans symétrie et des bords à courbure oscillante.
- **Termes d'ordre supérieur.** Au-delà de R^{−1/2}, la décomposition n'est pas étudiée ; le test sans paramètre libre (χ²/n = 1,16) n'en exige pas à la précision actuelle. D'après le calcul de Δ(R), le terme suivant vaut moins de 0,0015 sur C\_eff pour R ≥ 10³, et moins de 0,0003 pour R ≥ 10⁴.

**Questions ouvertes.**

1. Confirmer que le coefficient B de Huxley (2014) pour les sommets vaut (c₀/2π)·Ω(K), ce que le modèle local reproduit.
2. Chercher une forme close, ou au moins une expression plus compacte, pour la partie intermédiaire (a ≈ 0,3 à 3) de l'intégrale de C∞, à l'image de la constante c₀ de Balog et Deshouillers.
3. Démontrer l'équirépartition de la projection du vecteur partenaire, utilisée dans le calcul de la contribution des grands a à B.
4. Démontrer le passage de ℤ² au réseau aléatoire pour le déficit de périmètre (équirépartition des longs horocycles, avec un contrôle des excursions dans la pointe).
5. Étendre à la dimension 3 : déficit de surface de l'enveloppe entière, avec un invariant candidat de la forme ∫K^{a} dS, K la courbure de Gauss.

## Références

- A. Balog, I. Bárány, [On the convex hull of the integer points in a disc](https://elischolar.library.yale.edu/cowles-discussion-paper-series/1205/), Cowles Foundation DP 962, 1990.
- I. Bárány, D. G. Larman, [The convex hull of the integer points in a large ball](https://link.springer.com/article/10.1007/s002080050217), Math. Annalen, 1998.
- I. Bárány, [Extremal problems for convex lattice polytopes: a survey](https://users.renyi.hu/~barany/cikkek/psurveyAMS.pdf), Contemporary Mathematics.
- I. Bárány, M. Prodromou, [On maximal convex lattice polygons inscribed in a plane convex set](https://link.springer.com/article/10.1007/BF02773612), Israel J. Math., 2006.
- A. Balog, J.-M. Deshouillers, constante moyenne 3,453…, citée dans [The integer hull of the set {xy ≥ N}](https://arxiv.org/html/2602.06897), arXiv 2602.06897.
- M. N. Huxley, [The convex hull of the lattice points inside a curve](https://link.springer.com/article/10.1007/s10998-014-0024-5), Period. Math. Hungar. 68(1), 100–118, 2014.
- J.-O. Lachaud, D. Coeurjolly, T. Roussillon, [Geometry of Gauss digitized convex shapes](https://perso.liris.cnrs.fr/david.coeurjolly/publication/dgmm-25-geometryconvex/dgmm-25-geometryconvex.pdf), DGMM 2025.
- T. Roussillon et al., [Multigrid convergence of discrete geometric estimators](https://perso.liris.cnrs.fr/tristan.roussillon/static/file/2012bookchapter.pdf), 2012.
- F. de Vieilleville, J.-O. Lachaud, F. Feschet, [Convex digital polygons, maximal digital straight segments and convergence of discrete geometric estimators](https://link.springer.com/article/10.1007/s10851-007-0779-x), JMIV, 2007.
- S. Finch, [Convex lattice polygons](https://oeis.org/A249455/a249455.pdf), 2003.

## Annexe : reproductibilité

Toutes les valeurs de la §5 sont produites par les scripts Python ci-dessous (NumPy, SciPy, pandas). Ils ont été exécutés sur deux machines (Linux, et Windows avec Python 3.14 et NumPy 2.x), avec des résultats identiques au chiffre près, les graines aléatoires ne dépendant pas de la machine.

| Script | Rôle | Sortie |
| --- | --- | --- |
| balayage\_translations.py | 11 formes × 26 rayons (100 à 30 000) × 64 tirages ; invariants exacts | balayage\_translations.csv, invariants.csv |
| analyse\_translations.py | C\_eff par fenêtre avec erreurs, dispersion, rapport des sommets | C\_eff\_par\_R.csv |
| balayage\_grands\_R.py | 4 formes × 9 rayons (30 000 à 10⁶) × 128 tirages, calcul parallèle avec reprise | balayage\_grands\_R.csv |
| analyse\_grands\_R.py | comparaison des lois de correction : prédiction à l'aveugle, χ², pente locale | tableaux console |
| decomposition\_aretes.py | part de C\_eff par tranche de longueur d'arête (disque) | decomposition\_aretes.npy |
| queue\_aretes.py | queue des arêtes longues et borne ℓ'max = 2√2·R^{1/6} | queue\_aretes.npy |
| modele\_local.py | modèle local de Haar par Monte-Carlo : densité de sommets, C∞, queue en y^{−1/4} | modele\_local.npy, journal par tranche |
| formule\_integrale.py | formule intégrale de C∞ (§3.4), calcul autonome en une minute | journal console |
| queue\_cible.py, assemblage\_C.py | configurations à deux colonnes (vérification de L = 1/9) ; assemblage de référence de C∞ | journal de calcul |
| constante\_B.py | B des petits a par la fonction zêta d'Epstein, avec contrôle de la régularisation (mpmath requis) | tableau console |
| constante\_B\_grand\_a.py | B des grands a par la somme sur les vecteurs partenaires entiers | tableau console |
| test\_classes\_w0.py, delta\_R.py | vérification classe par classe sur le disque réel ; écart somme discrète − intégrale | tableaux console |
| identite\_pick.py | vérification de l'identité de Pick dans le modèle local : ⟨g⟩ = n\_∂/2 (§3.6) | tableau console |
| verif\_rapport.py | contrôles ponctuels (disque en translation, I₄/₃ exact, R = 30 000) | tableau console |

Contrôles intégrés : P\_hull ≤ P\_cont pour chaque calcul ; identité avec une enveloppe calculée sur grille pleine pour 33 cas (réseau décalé) ; pour le modèle local, identité avec une énumération brute sur 150 réseaux et indépendance vis-à-vis de la largeur de fenêtre.

**Erreurs rencontrées et corrigées.** Elles sont listées parce que chacune produisait un résultat plausible mais faux.

1. **Grille trop petite.** La forme était rognée sur le carré \[−R−1, R+1\]², ce qui créait un faux plateau du déficit relatif pour p > 2.
2. **Moyenne sur les seules rotations.** Pour le disque, les tirages étaient tous identiques, d'où un faux creux à R = 1 000 et une valeur trop basse.
3. **Seuil dans l'équation du second degré (modèle local).** Pour un vecteur court presque vertical, l'équation était traitée à tort comme linéaire, ce qui produisait de rares valeurs aberrantes ; détecté grâce à la réplication sur la seconde machine.
4. **Lois d'extrapolation.** Une correction en R^{−1/3}, puis en 1/ln R, avaient d'abord été retenues ; les données jusqu'à 10⁶ et le mécanisme de la §3.3 conduisent à R^{−1/6}.
5. **Calcul de B.** Une première prédiction (0,387) moyennait chaque direction d'arête sur toutes les orientations du bord, alors qu'une arête de direction w₀ n'apparaît qu'au point où la tangente lui est parallèle ; la seconde (0,2315) omettait la contribution des grands a. Le test classe par classe sur le disque réel a permis de localiser les deux erreurs.
