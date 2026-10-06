# DYN - Finalisation Zenodo v0.1b

Patrick Royer - chercheur indépendant, sans affiliation institutionnelle.

Accès public ; Licence Dynagénèse Éthique (DEL) personnalisée, version 1.1.

Les développements ont été réalisés avec l'aide d'IA. ChatGPT et Codex (OpenAI) participent au développement et aux contrôles documentés ; ChatGPT assure l'intégration technique et documentaire finale sous la direction de Patrick Royer. Contributions et outils gratuits : metadata/NOTICE_IA_ET_OUTILS.md.

Cette version distincte conserve les données brutes, les résultats et la portée scientifique de v0.1a. Elle renseigne l'auteur, ajoute l'attribution des IA et outils, adapte la DEL, précise la justification de l'angle réduit par la symétrie du réseau carré et corrige le renvoi C20 vers la fiche disponible. La loi à deux termes demeure conjecturale et P2 NOT_RUN.

## Fichiers à présenter

- manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.pdf et source Markdown.
- licenses/DYN-LICENCE_DEL_v1_1.pdf, texte français et traduction anglaise.
- metadata/DESCRIPTION.md, notice des contributions et champs de licence personnalisée.
- data/raw, code/reanalysis et evidence : noyau auditable et limites.

## Reproduction numérique depuis la racine

    python code/reanalysis/reanalyse_zenodo_candidate.py --root . --output evidence/reanalysis_reproduced

    python code/reanalysis/recompute_superellipse_stats.py 1.75 --data-dir data/raw

Ces scripts sont identiques à v0.1a. Aucun générateur n'est nécessaire à ces recalculs. Le script documentaire code/documentation/build_publication_pdfs.py régénère les deux PDF ; ses dépendances sont distinctes du calcul scientifique.

## Licence

Lire licenses/DEL_v1_1_FR.md et licenses/LICENSE_SCOPE.md. La licence publique comporte des restrictions d'usage et ne se présente pas comme une licence open source OSI. Les logiciels tiers restent sous leurs licences. Les anciennes propositions CC BY/MIT et le modèle DEL original sont conservés comme historique, sans être la notice active de cette version.

## Dépôt

Les champs humains confirmés et la licence personnalisée sont préparés dans metadata/. Le JSON est documentaire, pas une charge API validée. Le propriétaire doit rattacher le DOI réel et adapter la date de publication au dépôt. Cette préparation n'a effectué aucun envoi ou publication externe.
