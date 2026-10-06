# Corrective documentaire Zenodo v0.1a — déficit de périmètre

Corrective distincte préparée le 4 octobre 2026 à partir de la candidate locale v0.1. Elle n'a pas été publiée sur Zenodo. L'archive v0.1 d'origine reste inchangée et les fichiers de `data/raw/` sont recopiés à l'identique.

## Objet et portée

Cette version applique les corrections F01–F06 de la contre-vérification : commandes de reproduction, mesure de randomisation, espérance de l'estimateur de pi, rendu mathématique du PDF, statut exploratoire de la superellipse p=3 et portée du registre C08.

Aucun générateur d'enveloppes n'a été exécuté, aucune campagne P2 n'a été lancée et aucune publication externe n'a été effectuée. La loi à deux termes reste une conjecture numérique.

## Résultats reproduits

Sur les 92 moyennes groupées, le modèle fixé `C_inf=0.7192`, `B=0.2894` donne `chi2=106.843594`, `chi2/n=1.161343` et `max|z|=2.693532`. L'ajustement pondéré diagonal donne `C_inf=0.7204066 +/- 0.0005488`, `B=0.2960535 +/- 0.0032784`, `chi2=101.993063`, `ddl=90`, soit `chi2/ddl=1.133256`.

Les trois formes à courbure strictement positive constituent le noyau principal (69 moyennes). La superellipse p=3, dont la courbure s'annule aux quatre points axiaux, est conservée comme extension exploratoire (23 moyennes). Les erreurs indiquées sont nominales : `SE=SD(ddof=1)/sqrt(n)` sous hypothèse de covariance diagonale. Elles n'incluent ni covariance croisée, ni sélection de modèle, ni erreur systématique géométrique, ni incertitude totale des constantes.

## Reproduction depuis la racine

Prérequis : Python, NumPy, pandas et SciPy. Ouvrir un terminal dans la racine de cette corrective, puis exécuter exactement :

```powershell
python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced
```

Pour le contrôle descriptif p=1,75 :

```powershell
python code/reanalysis/recompute_superellipse_stats.py 1.75 --data-dir data/raw
```

Ces deux commandes lisent seulement les données existantes. Leurs contrôles d'exécution v0.1a sont consignés dans `evidence/reproduction_v0_1a/`.

## Organisation

- `manuscript/` : source Markdown v0.1a et PDF recomposé.
- `evidence/` : corrective F01–F06, registre des affirmations, limites, reçus de reproduction et réanalyse.
- `code/reanalysis/` : scripts légers de lecture, regroupement et ajustement.
- `code/original/` : scripts historiques conservés pour provenance ; ne pas les relancer par défaut.
- `code/corrected_geometry_v0_3/` : contrôles géométriques finis archivés, sans portée asymptotique.
- `data/raw/` : données sources conservées octet pour octet depuis v0.1.
- `metadata/` : brouillon sans identité ni affiliation inférée ; champs humains bloquants laissés vides.
- `SHA256SUMS.txt` : empreintes de tous les fichiers du paquet, sauf lui-même.

## Avant un dépôt public

Le déposant doit encore fournir une identité vérifiée, décider s'il souhaite indiquer une affiliation ou un ORCID, confirmer les droits et licences, et contrôler le statut du DOI dans son brouillon Zenodo connecté. Le statut public exact du DOI reste inconnu. Le fichier de métadonnées ne doit pas être envoyé tel quel tant que ces champs ne sont pas complétés par la personne autorisée.

