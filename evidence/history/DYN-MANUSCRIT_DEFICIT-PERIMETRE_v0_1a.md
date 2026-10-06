# Déficit de périmètre des enveloppes entières randomisées

## Conjecture à deux termes, résultats numériques et applications à la mesure des courbes convexes

**Corrective documentaire distincte v0.1a — 4 octobre 2026**  
**Auteur :** non renseigné dans cette corrective ; aucune identité, affiliation ou ORCID n'est inféré.  
**Statut :** note de recherche numérique. La loi à deux termes est une conjecture, non un théorème.

## Résumé

Soit $K$ un corps convexe du plan et soit $P_R$ l'enveloppe convexe des points d'un réseau entier aléatoirement tourné et translaté qui appartiennent à $RK$. On étudie le déficit moyen de périmètre

$$
\Delta P_K(R)=P(RK)-\mathbb{E}\!\left[P(P_R)\right].
$$

Un résultat publié de Ngoc et Reitzner fixe l'ordre dominant $R^{-1/3}$, via le déficit de largeur moyenne. La présente note propose, sans le démontrer, le raffinement

$$
\Delta P_K(R)
=C_\infty I_{4/3}(K)R^{-1/3}
-B I_{3/2}(K)R^{-1/2}
+o\!\left(R^{-1/2}\right),
$$

où

$$
I_\alpha(K)=\int_{\partial K}\kappa(s)^\alpha\,ds,
\qquad
C_\infty=0{,}7192,
\qquad
B=0{,}2894.
$$

La réanalyse des deux CSV historiques, sans nouvelle génération, retrouve sur 92 moyennes groupées : $\chi^2=106{,}8436$, $\chi^2/n=1{,}1613$ et $\max|z|=2{,}6935$ pour les constantes fixées. L'ajustement pondéré diagonal donne $C_\infty=0{,}7204066\pm0{,}0005488$ et $B=0{,}2960535\pm0{,}0032784$, avec $\chi^2/90=1{,}1333$. Les incertitudes sont nominales et ne couvrent ni la sélection de modèle, ni une covariance croisée, ni les erreurs systématiques.

Les trois formes à courbure strictement positive constituent le noyau principal de 69 moyennes. La superellipse $p=3$, dont la courbure s'annule aux points axiaux, est conservée comme extension exploratoire de 23 moyennes. L'agrégat de 92 points ne doit donc pas être présenté comme quatre exemples relevant tous de la même classe principale.

**Mots-clés :** enveloppe convexe entière ; réseau randomisé ; déficit de périmètre ; intégrales de courbure ; géométrie numérique ; conjecture.

<!-- pagebreak -->

## 1. Statut mathématique et apport

Trois niveaux de preuve doivent rester séparés.

**Résultat publié.** En dimension deux, le travail de Ngoc et Reitzner établit des bornes de même ordre $R^{-1/3}$ pour le déficit moyen de largeur sous ses propres hypothèses de régularité. Par la relation entre largeur moyenne et périmètre dans le plan, ce résultat fixe l'ordre dominant pertinent. Il ne fournit ni les valeurs $C_\infty$ et $B$ proposées ici, ni l'expansion à deux termes ci-dessus.

**Résultat numérique reproduit.** La chaîne « CSV existants, invariants, regroupement, modèle diagonal » reproduit les résultats centraux. Elle atteste un calcul sur un corpus fini ; elle n'est ni une preuve asymptotique, ni un test confirmatoire préenregistré.

**Conjecture.** L'universalité des coefficients, la forme exacte du second terme et la petitesse uniforme du reste ne sont pas démontrées. Les mécanismes locaux, intégraux et arithmétiques donnent une motivation calculable, pas un théorème global.

La classe prudente retenue pour la conjecture est celle des corps convexes plans à bord $C^2$, strictement convexes, avec courbure positive et bornée. Cette classe de travail ne doit pas être confondue avec la classe exacte du théorème d'ordre cité, qui possède ses propres hypothèses.

## 2. Mesure de randomisation

Une écriture sans ambiguïté consiste à fixer la forme $RK$ et à randomiser le réseau :

$$
L_{\rho,t}=\rho\!\left(\mathbb{Z}^2+t\right),
\qquad
t\sim\operatorname{Unif}\!\left([0,1)^2\right),
\qquad
\rho\sim\operatorname{Haar}\!\left(SO(2)\right),
$$

avec $t$ et $\rho$ indépendants. L'enveloppe est alors

$$
P_R(\rho,t)=\operatorname{conv}\!\left(RK\cap L_{\rho,t}\right).
$$

Si l'on préfère écrire $L_{\rho,t}=\rho\mathbb{Z}^2+u$, alors $u=\rho t$ est, conditionnellement à $\rho$, uniforme dans la cellule tournée $\rho[0,1)^2$. Il ne faut pas déclarer $u$ uniforme dans un carré fixe indépendant de $\rho$ dans cette écriture.

En appliquant la rotation inverse, on obtient l'écriture équivalente où le réseau reste fixe et la forme est randomisée :

$$
\widetilde P_R(g,v)
=\operatorname{conv}\!\left((gRK+v)\cap\mathbb{Z}^2\right),
\qquad
g=\rho^{-1},
\qquad
v=-t\pmod{\mathbb{Z}^2}.
$$

Ici $g$ est uniforme sur $SO(2)$ et $v$ est uniforme dans le carré fixe $[0,1)^2$, indépendamment. Ce transport inverse est la convention effectivement compatible avec les translations des scripts. La correction porte sur la notation de la mesure, pas sur les valeurs simulées.

Pour les formes à symétrie d'ordre quatre utilisées dans le noyau numérique, les scripts échantillonnent l'angle modulo cette symétrie, sur $[0,\pi/2)$. Cette réduction est adaptée à ces formes ; elle ne définit pas à elle seule la loi générale pour un corps arbitraire.

<!-- pagebreak -->

## 3. Homogénéité et conjecture à deux termes

Sous l'homothétie de rapport $R$, la courbure devient $\kappa/R$ et l'élément d'arc devient $R\,ds$. Ainsi,

$$
I_\alpha(RK)=R^{1-\alpha}I_\alpha(K).
$$

Cette identité explique les puissances dimensionnelles des deux termes proposés. La constante effective est

$$
C_{\mathrm{eff}}(K,R)
=\frac{\Delta P_K(R)R^{1/3}}{I_{4/3}(K)}.
$$

La conjecture s'écrit alors

$$
C_{\mathrm{eff}}(K,R)
=C_\infty
-B\,\frac{I_{3/2}(K)}{I_{4/3}(K)}R^{-1/6}
+o\!\left(R^{-1/6}\right).
$$

Pour le disque unité, $I_{4/3}=I_{3/2}=2\pi$, d'où

$$
C_{\mathrm{eff}}(R)
=C_\infty-BR^{-1/6}+o\!\left(R^{-1/6}\right).
$$

La convergence en $R^{-1/6}$ est lente. Un plateau observé à des rayons modérés ne suffit donc pas à identifier séparément le terme dominant et le terme correctif.

## 4. Origine et portée des constantes candidates

Le calcul historique de $C_\infty$ part d'un modèle local parabolique et assemble une intégrale numérique. Les tableaux modele_local.npy et formule_integrale.npy conservent des sorties historiques. Ils rendent la valeur voisine de $0{,}7192$ traçable, mais le transfert uniforme du modèle local au réseau randomisé global n'est pas démontré.

Le coefficient $B$ est historiquement décomposé en deux contributions numériques proches de $0{,}2315$ et $0{,}0579$. Leur somme rend $0{,}2894$ traçable. L'exhaustivité des classes de configurations et le contrôle asymptotique du reste ne sont pas établis.

La notation $C_\infty=0{,}7192\pm0{,}0003$ rencontrée dans les calculs historiques décrit une estimation numérique ou systématique, et non l'erreur-type d'un échantillon aléatoire. Aucune incertitude probabiliste totale n'est documentée pour $B$.

La présente réanalyse atteste seulement que $C_\infty=0{,}7192$ et $B=0{,}2894$ sont tenus fixes pendant le recalcul du $\chi^2$. Elle n'atteste pas qu'ils ont été élaborés historiquement sans consultation des mêmes données. Cette limitation borne explicitement l'affirmation C08.

<!-- pagebreak -->

## 5. Données et protocole

Le fichier balayage_translations.csv contient 18 304 lignes et le fichier balayage_grands_R.csv 4 608 lignes. La sélection centrale retient quatre formes présentes dans la campagne grands rayons et les groupes tels que $R\ge1230$. Le rayon $R=30\,000$, présent dans les deux campagnes, regroupe 64 et 128 tirages distincts.

Le corpus sélectionné comprend 8 448 tirages, regroupés en 92 moyennes :

| Statut | Forme | Points groupés | Classe |
|---|---:|---:|---|
| principal | disque | 23 | courbure strictement positive |
| principal | ellipse 2:1 | 23 | courbure strictement positive |
| principal | ellipse 5:1 | 23 | courbure strictement positive |
| exploratoire | superellipse $p=3$ | 23 | courbure nulle aux quatre axes |

Pour chaque groupe $j$, l'erreur-type nominale est

$$
\operatorname{SE}_j=\frac{s_j}{\sqrt{n_j}},
$$

où $s_j$ est l'écart-type empirique calculé avec $n_j-1$ au dénominateur. Les tailles sont $n_j=64$, $128$ ou $192$ au rayon commun. Les 92 observations analysées sont donc des moyennes de groupes, pas des tirages bruts.

Le calcul descriptif utilise

$$
\chi^2
=\sum_j
\left(
\frac{C_{\mathrm{eff},j}-C_{\mathrm{pred},j}}
{\operatorname{SE}_j}
\right)^2,
$$

avec une matrice de covariance diagonale. Les fichiers ne permettent pas d'estimer les covariances entre rayons ou entre formes. Les erreurs-types mesurent la dispersion Monte-Carlo dans chaque groupe ; elles n'intègrent ni erreur de modèle, ni sélection du seuil, ni incertitude géométrique systématique.

## 6. Résultats reproduits

| Modèle | Paramètres | $\chi^2$ | Dénominateur | Rapport |
|---|---:|---:|---:|---:|
| constantes fixées | $C_\infty=0{,}7192$ ; $B=0{,}2894$ | 106,843594 | $n=92$ | 1,161343 |
| ajustement libre diagonal | $0{,}7204066\pm0{,}0005488$ ; $0{,}2960535\pm0{,}0032784$ | 101,993063 | ddl $=90$ | 1,133256 |
| $C_\infty$ fixé | $B=0{,}2892607\pm0{,}0010969$ | 106,827473 | ddl $=91$ | 1,173928 |

Le résidu standardisé maximal du modèle fixé vaut $2{,}693532$. Dans l'ajustement libre, la corrélation nominale entre $C_\infty$ et $B$ vaut $0{,}942365$ : les deux coefficients sont difficilement séparables sur cette fenêtre d'échelles.

La décomposition descriptive du $\chi^2$ fixé est :

| Forme | $\chi^2$ | $\chi^2/n$ |
|---|---:|---:|
| disque | 40,1967 | 1,7477 |
| ellipse 2:1 | 27,0405 | 1,1757 |
| ellipse 5:1 | 20,1198 | 0,8748 |
| superellipse $p=3$, exploratoire | 19,4865 | 0,8472 |

Le noyau principal de 69 points donne donc $\chi^2=87{,}3571$ et $\chi^2/n=1{,}2660$. Cette stratification montre que l'accord agrégé n'est pas uniforme entre formes. Aucun test confirmatoire ni aucune conclusion asymptotique n'est tiré de cette décomposition.

## 7. Incertitudes et réserves statistiques

Les valeurs précédentes sont reproduites à partir des sorties disponibles, mais leur interprétation reste descriptive.

1. Le modèle, les coefficients candidats et le seuil $R\ge1230$ ont été développés dans le même programme de recherche que les données.
2. Les erreurs-types nominales supposent une covariance diagonale et ne couvrent pas la dépendance induite par les choix communs de code, d'invariants et de modèle.
3. Les erreurs systématiques de construction géométrique, de quadrature et de format flottant ne sont pas incluses dans les barres.
4. Plusieurs formes de correction ont été explorées ; l'incertitude de sélection de modèle n'est pas propagée.
5. La proximité d'un $\chi^2$ avec son dénominateur n'est pas une preuve de l'universalité asymptotique.

La contre-vérification autonome a contrôlé les tailles, le grain des campagnes, l'absence de valeurs manquantes ou non finies, les domaines des angles et translations, les déficits négatifs et les doublons. Les différences entre $P_{\mathrm{cont}}-P_{\mathrm{hull}}$ et la colonne de déficit restent au niveau des arrondis float64. Ses quadratures en double précision ne sont pas des certificats par intervalles.

## 8. Superellipse p=3 et autres extensions

Pour la superellipse $|x|^p+|y|^p\le1$ avec $p=3$, la courbure s'annule aux quatre points axiaux. Cette forme sort donc de la classe principale imposant une borne inférieure strictement positive sur la courbure. Ses 23 moyennes sont conservées pour continuité et comme test exploratoire de robustesse, jamais comme quatrième membre du noyau théorique homogène.

Pour $p=1{,}75$, le JSON brut contient 5 888 périmètres répartis en sept groupes : 1 024 tirages aux cinq premiers rayons, 512 à $R=300\,000$ et 256 à $R=1\,000\,000$. Le recalcul descriptif donne

$$
\chi^2_{\mathrm{sans\ couche}}=12{,}1043,
$$

et, pour l'exposant fixé $\gamma=16/33$,

$$
A=-0{,}06041\pm0{,}02366,
\qquad
\chi^2=5{,}5879
\quad
(\mathrm{ddl}=6).
$$

Avec $\gamma$ libre,

$$
\gamma=0{,}2802\pm0{,}2058,
\qquad
A=-0{,}01191\pm0{,}02197,
\qquad
\chi^2=4{,}7047
\quad
(\mathrm{ddl}=5).
$$

L'exposant est peu identifié. Ces chiffres dépendent légèrement du périmètre continu soustrait à grand rayon ; leurs dernières décimales ne sont pas certifiées. Le cas $p=1{,}5$ reste incomplet faute de JSON brut, et la campagne P2 porte le statut NOT_RUN.

<!-- pagebreak -->

## 9. Audit géométrique archivé

Le paquet geometry_audit_v0_3 documente des contrôles finis, dont un chemin rationnel limité et des tests de cohérence. Ces contrôles ne valident pas rétroactivement toutes les valeurs historiques du périmètre d'enveloppe et ne démontrent aucune expansion asymptotique.

Deux routines flottantes concordantes constituent un indice de robustesse, pas un certificat d'exactitude. Les scripts historiques sont conservés pour provenance ; ils ne doivent pas être relancés pour reproduire la présente corrective.

## 10. Applications conditionnelles

Les applications ci-dessous découlent algébriquement de la conjecture. Aucun gain pratique mesuré n'est revendiqué.

### 10.1 Périmètre d'une forme numérisée

Pour un pas de maille $h$ et une forme de taille caractéristique $a$, le grand paramètre est $R=a/h$. Une correction conditionnelle du biais moyen peut être construite à partir des deux termes. Sa pertinence dépend de la régularité de la forme, du modèle de randomisation et de la domination effective du reste aux résolutions utilisées.

### 10.2 Règle de résolution

Une règle de résolution fondée sur le seul terme dominant vise le biais moyen relatif. Elle ne borne ni l'erreur d'un tirage individuel, ni la variance, ni un risque de couverture. La présence du second terme et des erreurs systématiques impose de la présenter comme ordre de grandeur conditionnel.

### 10.3 Disque et estimation de pi

Considérons un disque de rayon physique $a$ observé sur une maille de pas $h$, et posons

$$
\pi_h=\frac{P_h}{2a},
\qquad
\varepsilon=\frac{h}{a}.
$$

Conditionnellement à la conjecture à deux termes et pour $a$ fixé,

$$
\pi-\mathbb{E}[\pi_h]
=\pi
\left[
C_\infty\varepsilon^{4/3}
-B\varepsilon^{3/2}
\right]
+o\!\left(\varepsilon^{3/2}\right).
$$

Une correction candidate est

$$
\widehat{\pi}_{\mathrm{corr}}
=\frac{\pi_h}
{1-C_\infty\varepsilon^{4/3}+B\varepsilon^{3/2}},
$$

à condition que le dénominateur soit non nul. Cette transformation vise l'annulation du biais moyen au niveau des termes modélisés. Elle conserve la dispersion entre tirages et ne garantit pas qu'une erreur individuelle soit de l'ordre du reste. Pour plusieurs tirages, la même formule peut être appliquée à la moyenne empirique de $\pi_h$, avec une incertitude Monte-Carlo estimée séparément.

## 11. Reproductibilité

Depuis la racine de la corrective v0.1a, la réanalyse centrale est exécutée par :

    python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced

Le contrôle descriptif $p=1{,}75$ est exécuté par :

    python code/reanalysis/recompute_superellipse_stats.py 1.75 --data-dir data/raw

Ces commandes remplacent respectivement l'option inexistante --candidate-root et la référence au fichier absent reanalyse_zenodo_main.py de la version v0.1. Leurs sorties et codes de retour sont archivés dans evidence/reproduction_v0_1a.

La réanalyse centrale lit seulement balayage_translations.csv et balayage_grands_R.csv, recalcule les intégrales $I_{4/3}$ et $I_{3/2}$, agrège les groupes et effectue les ajustements pondérés. Elle ne lance aucun générateur Monte-Carlo.

Les empreintes de data/raw dans v0.1a sont comparées à celles de la candidate v0.1. Le reçu de conservation doit indiquer zéro différence. SHA256SUMS.txt couvre toutes les pièces du paquet final sauf lui-même.

## 12. Corrections F01–F06

| Référence | Correction v0.1a | Statut |
|---|---|---|
| F01 | deux commandes remplacées par le point d'entrée existant et testées depuis la racine | appliquée |
| F02 | mesure du réseau tournée écrite avec $u=\rho t$ ; transport inverse explicité | appliquée |
| F03 | formule portée sur $\mathbb{E}[\pi_h]$ ; dispersion individuelle et dénominateur précisés | appliquée |
| F04 | égalité avec reste dans le résumé ; équations recomposées et contrôle visuel requis | appliquée |
| F05 | $p=3$ isolé comme extension exploratoire ; noyau principal de 69 points | appliquée |
| F06 | C08 borné aux constantes fixées pendant la réanalyse ; indépendance historique non attestée | appliquée |

## 13. Limites et domaine de validité

La conjecture n'est étayée numériquement que sur une fenêtre finie de formes et de rayons. Les singularités de courbure peuvent modifier les exposants ou introduire des contributions supplémentaires. La covariance croisée manque, le plan d'analyse n'a pas été préenregistré et aucun corpus confirmatoire tenu à l'écart n'est inclus.

Le DOI communiqué n'est pas présenté comme un enregistrement public vérifié. Son statut exact doit être contrôlé par le propriétaire du brouillon Zenodo. Les droits et licences doivent encore être arrêtés par la personne autorisée. La corrective ne renseigne aucune identité, affiliation ou ORCID par supposition.

## 14. Conclusion

Les sorties existantes corroborent la chaîne de réanalyse et montrent un accord descriptif avec la loi proposée sur le domaine observé. Elles ne démontrent ni les coefficients universels, ni le second terme, ni le reste. La formulation correcte est donc celle d'une conjecture numérique accompagnée d'un corpus auditable, d'incertitudes nominales et de réserves explicites.

## Références

1. B. H. Ngoc et M. Reitzner, « Expected mean width of the randomized integer convex hull », Mathematika 67 (2021), 422–433. DOI : 10.1112/mtk.12080.
2. I. Bárány et J. Matoušek, « The randomized integer convex hull », Discrete & Computational Geometry 33 (2005), 3–25. DOI : 10.1007/s00454-003-0836-1.
3. I. Bárány et D. G. Larman, « The convex hull of the integer points in a large ball », Mathematische Annalen 312 (1998), 167–181. DOI : 10.1007/s002080050217.

## Déclaration de disponibilité

Les données, scripts historiques, scripts de réanalyse, résultats, registres et empreintes sont inclus dans la corrective locale associée. Aucune publication Zenodo n'a été effectuée pendant sa préparation.
