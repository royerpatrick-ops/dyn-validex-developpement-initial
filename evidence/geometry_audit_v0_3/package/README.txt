DYN-VALIDEX — CORRECTIVE INTEGER HULL v0.3
Candidate géométrique distincte, 2 octobre 2026.

Cette version traite le défaut de conservation des points frontières et le
caractère incomplet du premier patch. Elle n'est pas une validation indépendante
déjà acquise, une certification de C, ni une démonstration de la conjecture.

À lire : NOTE_GEOMETRIE_v0_3.txt, VALIDATION_LOCALE.txt, PROMPT_AUDIT_CODEX.txt.
La note scientifique v0.2 est conservée dans history/ : son statut conjectural
et ses limites numériques restent inchangés. Les références à ses résultats
locaux concernent le paquet v0.2 conservé, pas des nouveaux calculs en v0.3.

Contenu
  candidate/ : les modules actifs v0.3, les tests et le contrôle géométrique.
  originals/ et document_original.pdf : les sources originales conservées.
  history/ : les paquets v0.2 et revue du patch, conservés sans modification.
  local_validation/ : résultats effectivement exécutés sur cette candidate.
  originals_sha256.json : empreintes des 89 fichiers sources originaux.
  candidate_changes_v0_2_to_v0_3.diff : changements depuis la candidate v0.2.
  package_manifest.json : empreintes de la livraison, sauf ce manifeste.

Principaux changements
  - geometry_columns.column_indices renvoie des indices entiers du réseau.
    Pour le disque (ellipse k=1 ou superellipse p=2), toutes les décisions de
    colonnes et d'appartenance utilisent des entiers exacts, même lorsque R/u/v
    ont été fournis sous forme de nombres float64.
  - Les autres formes utilisent une boîte englobante à bornes rationnelles,
    une recherche du minimum discret par convexité, puis des recherches des
    extrémités par appartenance filtrée. Aucun discriminant approché, aucune
    racine flottante ou minimisation dorée ne peut écarter une colonne valide.
  - Les décisions proches de la frontière sont reprises exactement lorsque
    possible, ou par intervalles à précision croissante. Un cas non résolu
    lève GeometryIndeterminate : aucun résultat partiel n'est renvoyé.
  - L'enveloppe utilise des orientations en entiers Python, sans débordement
    d'int64. Les longueurs restent évaluées en virgule flottante avec math.fsum.
  - Les commandes de campagne, d'estimation de C et de Monte Carlo par classes
    sont suspendues dans v0.3. Les noyaux mathématiques sont conservés pour examen.

API et sens des paramètres
  R, prm, theta, u, v acceptent int, float ou fractions.Fraction.
  Un float désigne sa valeur rationnelle binaire exacte. Pour demander le
  rationnel 1/3, utiliser Fraction(1,3), pas le float 1/3.
  Pour les autres formes que le disque, la rotation est celle de l'angle exact
  fourni ; ses sin/cos sont encadrés, et non simplement repris de math.sin/cos.
  Les translations doivent appartenir à [0,1) ; R>0 ; ellipse k>0 ; superellipse p>1.

  from fractions import Fraction
  from geometry_columns import column_indices, integer_hull_perimeter
  points = column_indices('ellipse', 1, 31, 0.314159, Fraction(1,3), Fraction(2,3))
  perimeter, vertices = integer_hull_perimeter(points)

  columns() reste disponible pour produire des coordonnées float translatées,
  destinées à la compatibilité/visualisation. Pour construire l'enveloppe active,
  utiliser column_indices(). hull_perim() refuse désormais les coordonnées float.
  Les conventions de périmètre exigent une enveloppe non dégénérée ; une enveloppe
  vide ou de dimension inférieure est refusée, sans valeur fabriquée.

Limites explicites
  - Par défaut, au plus 2 000 005 colonnes et des indices dans [-2^40,2^40].
    Un dépassement est refusé avant calcul ; max_columns est paramétrable.
  - Les intervalles sont raffinés à 50,100,200,400 chiffres décimaux par défaut.
    L'absence de décision finale provoque un arrêt explicite, pas une tolérance.
  - Le filtre suppose les opérations IEEE de NumPy et les fonctions d'intervalles
    de mpmath correctes. Les fonctions élémentaires de mpmath.iv sont documentées
    comme expérimentales : l'audit doit examiner cette dépendance. Les contrôles
    locaux ne constituent pas une preuve formelle de ce logiciel tiers.
  - Les cas exacts de disque sont indépendants de ce filtre et de la trigonométrie.
  - Cette approche privilégie la conservation des points et peut être plus lente
    que v0.2. Aucune nouvelle campagne complète n'a été exécutée.

Exécution de l'audit local, depuis la racine décompressée :

  python -m unittest discover -s candidate -p "test_*.py" -v
  python candidate/run_geometry_checks.py --out audit_codex_v03/geometry

Le dossier de sortie doit être neuf. Ne pas remplacer les résultats locaux.
Transmettre PROMPT_AUDIT_CODEX.txt à Codex avec le dossier complet. L'audit demandé
est géométrique et régressif ; aucun nouveau calcul de C n'est autorisé par ce prompt.
