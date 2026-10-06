# Dictionnaire de données

## CSV centraux

`raw/balayage_translations.csv` et `raw/balayage_grands_R.csv` partagent les colonnes :

- `forme` : famille (`ellipse` ou `super`).
- `param` : rapport d'axes pour l'ellipse; exposant `p` pour la superellipse.
- `R` : facteur d'échelle en unités du réseau.
- `echantillon` : identifiant du tirage dans le groupe.
- `theta` : rotation en radians.
- `u`, `v` : translation du réseau dans `[0,1)^2`.
- `P_cont` : périmètre continu de la forme dilatée.
- `P_hull` : périmètre calculé de l'enveloppe entière.
- `deficit_abs` : `P_cont - P_hull`.
- `m_sommets` : nombre de sommets de l'enveloppe.

`raw/invariants.csv` contient `P`, `Omega=I_(1/3)` et `I43=I_(4/3)` pour les onze formes de la campagne courte. `raw/C_eff_par_R.csv` est une agrégation historique, non la source primaire de la réanalyse indépendante.

Dans la convention des CSV, la forme est tournée par rapport au réseau fixe et `u,v` sont uniformes dans le carré fixe `[0,1)^2`. La description équivalente à forme fixe utilise un réseau `rho(Z^2+t)`; si elle est écrite `rho Z^2+u`, la translation `u=rho t` appartient conditionnellement à la cellule tournée.

## Tailles et incertitudes

Dans la réanalyse des 92 points, `n` est le nombre de tirages bruts par groupe forme-rayon : 64 dans la campagne courte, 128 dans la campagne grand R, et 192 après regroupement au rayon commun 30000. `C_eff_se` est `SD(ddof=1)/sqrt(n)`. La covariance croisée n'est pas disponible et est supposée nulle dans le calcul diagonal.

Le noyau à courbure strictement positive contient 69 moyennes (disque et ellipses 2:1 et 5:1). Les 23 moyennes de la superellipse p=3 sont conservées comme extension exploratoire, car sa courbure s'annule aux points axiaux.

## JSON p=1.75

`raw/superellipse_p1.75_singuliere.json` est un dictionnaire : clé externe = rayon sous forme de chaîne; valeur = dictionnaire `identifiant_echantillon -> P_hull`. Il y a sept rayons : 1024 valeurs pour chacun des cinq premiers, 512 à `R=300000` et 256 à `R=1000000`. Le périmètre continu et les invariants sont recalculés par `code/reanalysis/enveloppe_commune.py`.

## p=1.5

`raw/test_superellipse_p1.5_2026-10-03.md` est une synthèse agrégée, pas un jeu brut. Aucun JSON p=1.5 n'est présent. Les chiffres qu'il contient sont `PRESENT_UNVERIFIED` au niveau échantillon.

## NPY

- `raw/decomposition_aretes.npy` : dictionnaire objet, clés `10000`, `100000`, `1000000`; valeurs tuples de résumés d'arêtes.
- `raw/formule_integrale.npy` : tableau `float64`, forme `(29,4)`.
- `raw/modele_local.npy` : dictionnaire objet avec `tranches` `(25,7)`, et vecteurs `w`, `c`, `p` de longueur `172155`.
- `raw/queue_aretes.npy` : dictionnaire objet avec `bins` de longueur 14 et quatre histogrammes de longueur 13.

Les NPY objet nécessitent `numpy.load(..., allow_pickle=True)` et ne doivent être ouverts qu'à partir de cette archive de confiance.

## Valeurs dérivées

Les fichiers sous `evidence/reanalysis/` sont dérivés. `reanalysis_92_points.csv` contient `C_eff`, son écart-type, son erreur-type, les invariants, la prédiction fixée et le résidu standardisé. `reanalysis_summary.json` fixe les conventions et les matrices de covariance nominales.

## Données absentes

- JSON brut p=1.5 : absent.
- NPZ de localisation : absent.
- P2 : `NOT_RUN`, aucune donnée.
