# Corrections cumulées

## Corrections scientifiques et documentaires antérieures

1. L'ordre dominant R^(-1/3) publié par Ngoc et Reitzner est reconnu.
2. La loi à deux termes est qualifiée de conjecture soutenue numériquement.
3. C_inf=0.7192 et B=0.2894 sont des valeurs candidates.
4. Le modèle local n'est pas présenté comme une dérivation globale rigoureuse.
5. La concordance de routines géométriques sur des cas finis n'est pas une preuve générale.
6. Le sens des erreurs-types, de n et du modèle diagonal est explicité.
7. La reproduction sur une seconde machine reste non vérifiée.
8. Le cas p=1.5 reste incomplet faute de JSON brut.
9. La localisation reste non testée faute de NPZ.
10. Le retrait d'un ancien argument d'anisotropie ne prouve pas l'isotropie.
11. Pour la lemniscate, la cible d'une enveloppe est le périmètre de l'enveloppe convexe.
12. Les applications sont conditionnelles et sans gain mesuré revendiqué.
13. La nouveauté est bornée à une combinaison explicite et testable.
14. La campagne P2 reste NOT_RUN.

## Corrective documentaire v0.1a

15. **F01 — commandes.** Les deux appels erronés sont remplacés par la commande centrale : python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced. Le contrôle p=1.75 utilise explicitement --data-dir data/raw.
16. **F02 — randomisation.** Le réseau est écrit rho(Z^2+t), avec t uniforme dans la cellule fixe et rho distribuée selon Haar, indépendamment. Dans l'écriture rho Z^2+u, la translation conditionnelle appartient à la cellule tournée. Le transport inverse vers la forme randomisée est explicité.
17. **F03 — pi.** Le développement porte sur E[pi_h], non sur chaque tirage. L'estimateur corrigé requiert un dénominateur non nul et ne supprime pas la dispersion individuelle.
18. **F04 — rendu.** Le résumé utilise une égalité avec reste. Les équations affichées sont composées lors de la génération du PDF, puis les pages sont rendues et inspectées.
19. **F05 — p=3.** Les 69 points des trois formes à courbure strictement positive forment le noyau ; les 23 points p=3 sont exploratoires.
20. **F06 — C08.** Seul le maintien des constantes fixées pendant la réanalyse est vérifié. Leur indépendance historique à l'égard des mêmes données n'est pas attestée.
