# DYN-VALIDEX — développement initial

Archive scientifique associée au DOI [10.5281/zenodo.23110956](https://doi.org/10.5281/zenodo.23110956).

Auteur humain et responsable de la publication : **Patrick Royer**, chercheur indépendant sans affiliation institutionnelle.

## Portée scientifique exacte

Ce dépôt documente une **conjecture numérique à deux termes** pour le déficit moyen de périmètre de l’enveloppe convexe des points d’un réseau entier randomisé contenus dans un grand corps convexe plan. La contribution proposée concerne des coefficients de courbure, un terme correctif en `R^(-1/2)`, deux constantes candidates et leur évaluation numérique.

Le dossier contient une réanalyse, sans nouvelle simulation, de 92 moyennes archivées. Les 69 moyennes du disque et de deux ellipses constituent le noyau à courbure strictement positive ; les 23 moyennes de la superellipse `p=3` sont exploratoires.

La loi exacte à deux termes, l’universalité des coefficients et le reste annoncé **ne sont pas démontrés**. Les constantes restent candidates, l’audit géométrique a une portée finie, les incertitudes statistiques sont nominales sous covariance diagonale et la campagne P2 est au statut `NOT_RUN`. Les limites et pièces manquantes sont détaillées dans [`evidence/LIMITS_AND_MISSING.md`](evidence/LIMITS_AND_MISSING.md).

## Contenu principal

- [`manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.pdf`](manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.pdf) — rapport principal ;
- [`manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.md`](manuscript/DYN-MANUSCRIT_DEFICIT-PERIMETRE_v0_1b.md) — source éditable du rapport ;
- [`code/`](code/) et [`data/`](data/) — sources et données disponibles ;
- [`evidence/`](evidence/) — registres d’affirmations, contrôles, rapports et empreintes ;
- [`licenses/DYN-LICENCE_DEL_v1_1.pdf`](licenses/DYN-LICENCE_DEL_v1_1.pdf) — licence DEL 1.1, avec textes source dans [`licenses/`](licenses/) ;
- [`metadata/NOTICE_IA_ET_OUTILS.md`](metadata/NOTICE_IA_ET_OUTILS.md) — déclaration détaillée de l’assistance IA et des outils.

## Assistance par intelligence artificielle

Les développements ont été réalisés avec une assistance substantielle d’IA. ChatGPT et Codex (OpenAI) ont participé à l’exploration des formulations, au développement et à la correction de scripts, à l’analyse de données disponibles, à la rédaction, à des vérifications numériques et documentaires et à l’intégration du dossier, sous la direction de Patrick Royer. Les contrôles assistés par IA sont des contrôles bornés et documentés ; ils ne constituent ni une preuve générale ni une revue humaine extérieure.

## Statut de la relecture

**Aucune revue scientifique humaine extérieure n’est attestée dans ce dépôt.** Les audits ou contrelectures identifiés comme provenant de ChatGPT/Codex sont des revues par IA. Une relecture humaine indépendante est explicitement sollicitée dans [`REVIEW_REQUEST.md`](REVIEW_REQUEST.md).

## Reproductibilité et intégrité

Les scripts de réanalyse et de reproduction archivés sont fournis pour audit. La préparation de ce dépôt GitHub n’a relancé aucun calcul scientifique, aucune campagne P2 et aucun flux RNG. Les certificats, résultats et pièces gelées ont été copiés sans modification ; seuls ce README et la demande de relecture ont été ajoutés pour l’accueil GitHub.

## Licence

Les contenus couverts sont publiés sous la **Licence Dynagénèse Éthique (DEL) 1.1**, licence personnalisée comportant des restrictions d’usage et ne se présentant pas comme une licence open source approuvée par l’OSI. Les faits et concepts mathématiques, les composants tiers et les autorisations antérieures conservent leurs régimes propres. Lire [`licenses/LICENSE_SCOPE.md`](licenses/LICENSE_SCOPE.md).
