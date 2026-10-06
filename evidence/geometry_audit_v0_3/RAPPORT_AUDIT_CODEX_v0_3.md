# Audit indépendant DYN-VALIDEX — candidate INTEGER HULL v0.3

Date : 2026-10-02  
Paquet audité : `DYN-CORRECTIVE_INTEGER-HULL_v0_3`  
Périmètre : géométrie entière et régressions finies exclusivement.

## Verdict

**PASS pour G01–G16 dans le périmètre fini annoncé.** Aucun défaut bloquant de la corrective géométrique v0.3 n'a été reproduit. Les 33 tests prescrits, les contrôles structurés 1494/1494 et 780/780, ainsi que le harnais indépendant passent.

La voie exacte du disque corrige le défaut d'arrondi de v0.2 : aucun sommet manquant ou substitué n'est observé parmi les rayons, translations, frontières et grilles exhaustives contrôlés. Elle ne dépend ni de la trigonométrie flottante ni de `mpmath.iv`.

Pour les ellipses non circulaires et superellipses, le résultat est **validé sur le domaine fini testé**, mais ne constitue pas une certification globale : le dernier niveau de décision dépend de `mpmath.iv`, dont la documentation qualifie encore le support des intervalles d'expérimental.

Aucun calcul de constante `C`, aucune campagne Monte Carlo et aucun profil long n'ont été exécutés. Il n'est donc tiré aucune conclusion nouvelle sur `C`, le second terme, le reste ou la loi asymptotique globale.

## Environnement relevé

- Exécutable : `C:\Users\dynagenese\Desktop\pi utile et efficace\DYN-CORRECTIVE_INTEGER-HULL_v0_2\DYN-CORRECTIVE_INTEGER-HULL_v0_2\.venv_audit\Scripts\python.exe`
- CPython 3.14.5, Windows 11 `10.0.26300`
- NumPy 2.3.5, SciPy 1.17.0, mpmath 1.4.1
- Différence par rapport à la validation soumise : Python 3.14.5/Windows au lieu de Python 3.12.14/Linux ; versions des trois bibliothèques identiques.
- `PYTHONDONTWRITEBYTECODE=1` a été imposé aux exécutions.

Le détail machine est conservé dans `environment.json` et les commandes exactes dans `commands.txt`.

## Reproductions prescrites

1. `python -m unittest discover -s candidate -p "test_*.py" -v`
   - code retour 0 ; 33 tests exécutés ; `OK` ; 7,109 s rapportées par `unittest`.
2. `python candidate/run_geometry_checks.py --out audit_codex_v03/geometry`
   - code retour 0 ; statut `PASS_LOCAL_GEOMETRY_REGRESSIONS` ;
   - rayons 2 à 250, six angles : 1494/1494 ;
   - 25 anciens rayons défaillants plus `R=1000`, cinq translations, six angles : 780/780 ;
   - 25 rejeux historiques : même nombre de sommets ; écart absolu maximal de périmètre `1.8189894035458565e-12`.

Les sorties standard, erreurs standard, codes retour et durées sont conservés séparément. Le faible écart de périmètre est cohérent avec le calcul final des longueurs en flottants ; il ne constitue pas une borne d'erreur certifiée.

## Contrôles indépendants

Le script `independent_checks.py` refuse d'écraser son résultat et emploie des références distinctes du parcours de colonnes de production : énumération complète de grilles, inégalités rationnelles directes, enveloppe de Jarvis et évaluations ponctuelles à 220–260 chiffres.

- Intégrité : 122 entrées du manifeste, 89 originaux, 12 modules candidats ; toutes les empreintes concordent.
- Historique : 134 fichiers de l'archive v0.2 et 18 fichiers de la revue v0.1 comparés octet par octet à leurs paquets conservés.
- Disques rationnels : 30 cas exhaustifs, 1 236 points intérieurs classés ; dénominateurs 3, 5, 7, 9, 11, 97 et 101.
- Flottants de disque : 18 cas autour d'ULP, translations binaires exactes et plus petit subnormal positif ; trois paires frontière/extérieur immédiat.
- Formes générales : 11 familles, trois configurations chacune, soit 33 cas ; comparaison de toutes les extrémités de colonnes et des enveloppes sur 1 628 points intérieurs.
- Frontières générales : 18 classifications avec rayons flottants adjacents autour d'un rayon critique calculé à 240 chiffres.
- Encadrements : 1 728 intervalles float64 et 400 intervalles `mpmath.iv` confrontés à une évaluation à 260 chiffres ; aucune violation observée.
- Indétermination : un cas naturellement ambigu à 15 chiffres lève `GeometryIndeterminate`, sans inclusion ni exclusion silencieuse.
- Ressources et enveloppe : budgets de colonnes et de coordonnées refusés, orientation testée jusqu'à `±2^40`, coordonnées float non vides refusées, périmètre dégénéré refusé.
- Sorties : une sortie géométrique existante est refusée et reste inchangée ; les trois CLI suspendues sortent avec le code 2 sans créer de résultat ; les imports sont sans effet de bord.
- Régressions non géométriques : quatre modules sont identiques à v0.2 ; `log_nodes`, `deterministic_tail` et `estimate` ont un AST identique ; aucun calcul de `C` n'est effectué.

## Mécanisme de correction

### Disque exact

`disk_column_indices` convertit `R`, `u` et `v` en rationnels exacts, prend un dénominateur commun `D`, puis transforme l'appartenance en une inégalité entière. Pour chaque colonne `kx`, le résidu

`(R D)^2 - (D kx - U)^2`

est entier. Sa racine entière inférieure est calculée par `math.isqrt`; les bornes basse et haute sont ensuite des divisions euclidiennes exactes. Aucun `sqrt` flottant, aucune troncature approchée et aucune tolérance ne peuvent déplacer la frontière. Une rotation est ignorée uniquement sur cette voie parce que le disque est invariant par rotation.

L'erreur v0.2 venait au contraire d'une borne continue calculée en flottant puis arrondie avant la décision d'appartenance. Au voisinage d'un entier, une erreur d'un ULP pouvait diminuer ou augmenter une extrémité de colonne ; l'enveloppe remplaçait alors un sommet attendu par un voisin. Dans v0.3, les ensembles `missing` et `extra` sont vides sur tous les contrôles prescrits.

### Formes générales

Le filtre float64 propage des intervalles élargis vers l'extérieur. Il accepte seulement si la borne supérieure est `<= 1` et refuse seulement si la borne inférieure est `> 1`. Les cas non décidés passent à une évaluation rationnelle exacte lorsqu'elle existe, sinon à `mpmath.iv` aux précisions 50, 100, 200 et 400. Une ambiguïté persistante provoque `GeometryIndeterminate` et interrompt l'appel complet.

La recherche par colonne trouve d'abord le minimum discret d'une fonction strictement convexe, puis les deux frontières par recherches monotones. Les extrémités rendues sont revérifiées et leurs voisins immédiats doivent être extérieurs. Enfin, l'enveloppe utilise des orientations en entiers Python ; la conversion flottante n'intervient que pour les longueurs.

## G01–G16

| Critère | Statut | Preuve principale |
|---|---|---|
| G01 | PASS | Manifestes, 89 originaux, 12 candidats et deux archives historiques concordants. |
| G02 | PASS | Sémantique rationnelle exacte des int/float/Fraction ; 10 entrées invalides refusées. |
| G03 | PASS | Bornes du disque conformes sur 30 grilles rationnelles et 18 cas flottants. |
| G04 | PASS | Frontières incluses, extérieurs immédiats exclus ; 18 cas génériques proches de la frontière conformes. |
| G05 | PASS | 1494/1494 pour les rayons 2–250 à six angles. |
| G06 | PASS | 780/780 pour les 25 anciens rayons, `R=1000`, cinq translations et six angles. |
| G07 | PASS | Les anciens cas rationnels et `R=39,41` passent ; aucun sommet manquant ou supplémentaire. |
| G08 | PASS | Énumération rationnelle complète et enveloppe de Jarvis sur 30 cas. |
| G09 | PASS | Onze familles, rotations et translations : 33/33 cas indépendants. |
| G10 | PASS | Toutes les extrémités de colonnes comparées, pas seulement les sommets finaux. |
| G11 | PASS | 2 128 encadrements confrontés à 260 chiffres, sans violation ; portée empirique seulement. |
| G12 | PASS | Reprises exactes et par intervalles exercées ; indétermination naturelle arrêtée. |
| G13 | PASS | Symétrie, subnormal, budgets de colonnes et de coordonnées contrôlés. |
| G14 | PASS | Orientations entières à grande amplitude ; coordonnées float non vides refusées. |
| G15 | PASS | 33 tests ; noyaux invariants/K/B et fonctions numériques non géométriques préservés. |
| G16 | PASS | Sorties neuves, refus d'écrasement, imports propres, campagnes et `C` suspendus. |

## Réserve `mpmath.iv`

La documentation de contexte de mpmath décrit l'arithmétique d'intervalles comme assurant en principe le suivi rigoureux des erreurs, mais précise que ce support reste expérimental et que de nombreuses fonctions ne prennent pas encore correctement les intervalles en charge. La documentation trigonométrique donne explicitement des exemples de `iv.sin` et `iv.cos`.

Sources consultées :

- <https://mpmath.org/doc/current/contexts.html>
- <https://mpmath.org/doc/current/functions/trigonometric.html>

Les 400 contrôles directs de `mpmath.iv` et les 1 728 contrôles du filtre n'ont révélé aucune non-inclusion dans les opérations utilisées ici. Cela renforce la confiance pour la version et le domaine testés, sans transformer cette dépendance expérimentale en preuve formelle. Une correction frontière pleinement certifiée pour toutes les formes demanderait soit un backend d'intervalles audité/formellement garanti pour `sin`, `cos` et les puissances, soit des bornes mathématiques indépendantes accompagnées d'une preuve d'inclusion.

## Observations non bloquantes

1. `max_columns=True` est accepté comme entier Python et impose en pratique un budget de 1. Il n'autorise aucun dépassement, mais un rejet explicite des booléens rendrait le contrat de type homogène avec les paramètres géométriques.
2. `integer_hull_vertices` retourne immédiatement une enveloppe vide pour un tableau float vide avant de vérifier son type. Les tableaux float non vides sont bien refusés ; aucune coordonnée flottante n'est utilisée dans ce cas vide.

Ces deux points sont de robustesse d'API, sans contre-exemple géométrique ni contournement de ressource observé. Aucun patch séparé n'est proposé pour cette soumission.

## Conclusions séparées

- **Disque rationnel exact : PASS.** La correction de frontière est fiable dans le modèle d'entrée annoncé et supprime le mécanisme d'arrondi responsable des sommets manquants/substitués de v0.2.
- **Autres ellipses et superellipses : PASS sur les tests finis, certification globale INCONCLUSIVE.** Aucun contre-exemple n'est trouvé ; la garantie finale reste conditionnelle à `mpmath.iv` et aux hypothèses IEEE documentées.
- **Constante `C` : INCONCLUSIVE / non évaluée.** Aucun calcul nouveau n'a été autorisé ou effectué.
- **Théorie asymptotique globale : INCONCLUSIVE / non évaluée.** Transport vers Haar, second terme et reste ne sont ni prouvés ni infirmés par cette corrective.

## Recommandation

La v0.3 peut être retenue comme corrective géométrique distincte de v0.2 pour la voie exacte du disque et comme amélioration prudente, à arrêt explicite, pour les autres formes. Elle ne doit pas être présentée comme une certification de `C`, une preuve de la loi asymptotique ou une validation formelle de `mpmath.iv`.
