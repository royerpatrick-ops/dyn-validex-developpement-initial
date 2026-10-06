# Fiche corrective documentaire v0.1a

Date : 4 octobre 2026

## Base et invariants de conservation

- Base : DYN-CANDIDATE_ZENODO_DEFICIT-PERIMETRE_v0_1.zip
- SHA-256 de la base : d95c54ff82366f910ecfdd2df949c27f42b8a4201df70e7cfeaea529f2451407
- Contre-vérification : DYN-CONTRE-VERIFICATION_ZENODO_DEFICIT-PERIMETRE_v0_1.zip
- La candidate v0.1 et son ZIP sont conservés.
- Les fichiers de data/raw sont recopiés octet pour octet ; la comparaison SHA-256 finale est consignée dans DATA_PRESERVATION_CHECK.json.
- Aucun générateur d'enveloppes, aucune campagne P2 et aucune publication externe.

## F01 — commandes

Les deux appels défectueux de v0.1 sont remplacés par :

    python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced

Le contrôle p=1.75 s'exécute par :

    python code/reanalysis/recompute_superellipse_stats.py 1.75 --data-dir data/raw

Les deux commandes sont exécutées depuis la racine v0.1a. Les codes de retour, sorties et erreurs sont archivés sous evidence/reproduction_v0_1a.

## F02 — mesure de randomisation

À forme fixe, le réseau est L=rho(Z^2+t), avec t uniforme dans [0,1)^2 et rho Haar sur SO(2), indépendamment. Si L est écrit rho Z^2+u, alors u=rho t est conditionnellement uniforme dans rho[0,1)^2. Après transport inverse, la forme est tournée et translatée par rapport au réseau fixe ; la translation se trouve alors dans le carré fixe.

## F03 — espérance de pi_h

Le développement porte sur pi-E[pi_h]. L'estimateur corrigé est défini seulement pour un dénominateur non nul. Il vise le biais moyen et ne garantit aucune erreur individuelle de l'ordre du reste.

## F04 — PDF

Le résumé emploie une égalité avec reste o(R^-1/2), non une équivalence trompeuse à deux termes. Les formules affichées sont composées et le PDF fait l'objet d'un rendu page par page et d'un contrôle de glyphes.

## F05 — p=3

Les 69 moyennes du disque et des ellipses 2:1 et 5:1 forment le noyau à courbure strictement positive. Les 23 moyennes de la superellipse p=3 sont conservées comme extension exploratoire, puisque sa courbure s'annule aux points axiaux.

## F06 — C08

C08 atteste uniquement que C_inf=0.7192 et B=0.2894 sont maintenus fixes pendant le recalcul du chi2. Leur indépendance historique par rapport aux mêmes données n'est pas démontrée.

## Résultats statistiques reproduits

- Modèle fixé : chi2=106.843594433239 ; n=92 ; chi2/n=1.161343417753 ; max|z|=2.693532164134.
- Ajustement diagonal libre : C_inf=0.720406578490 +/- 0.000548761799 ; B=0.296053488200 +/- 0.003278351881 ; chi2=101.993063469903 ; ddl=90.
- Corrélation nominale des paramètres : 0.942364952306.
- C_inf fixé : B=0.289260726176 +/- 0.001096893118 ; chi2=106.827472750235 ; ddl=91.
- Noyau principal : n=69 ; chi2 fixé=87.357107987253.
- Extension p=3 : n=23 ; chi2 fixé=19.486486445986.

Les signes +/- désignent des erreurs-types nominales sous covariance diagonale. Ils n'incluent ni covariance croisée, ni choix de modèle ou de seuil, ni incertitude géométrique systématique, ni incertitude totale des constantes.

## Métadonnées

Le tableau creators est laissé vide. Aucune identité, affiliation ou ORCID n'est inféré. Le brouillon est explicitement bloqué pour publication tant que le déposant autorisé n'a pas renseigné le créateur, confirmé les droits et choisi les licences.
