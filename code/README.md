# Code

original/ contient les scripts historiques copiés sans modification. Les générateurs peuvent être coûteux et ne sont pas nécessaires à la finalisation documentaire v0.1b.

reanalysis/ contient les recalculs à partir des sorties existantes. Depuis la racine de la corrective :

    python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced

Contrôle descriptif p=1.75 :

    python code/reanalysis/recompute_superellipse_stats.py 1.75 --data-dir data/raw

corrected_geometry_v0_3/ contient les modules correctifs et tests de géométrie archivés. Leur audit a une portée finie et ne valide pas l'asymptotique.

requirements.txt décrit l'environnement de contrôle disponible en octobre 2026 ; il ne reconstitue pas nécessairement les environnements historiques complets.


documentation/ contient la composition des PDF v0.1b, développée avec l'assistance de ChatGPT/Codex. Ce rendu ne lance aucun calcul scientifique. Les dépendances de ce rendu figurent dans documentation/requirements-tested.txt.
