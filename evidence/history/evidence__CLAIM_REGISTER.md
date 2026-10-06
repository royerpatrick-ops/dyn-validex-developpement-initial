# Registre des affirmations

Le registre canonique est `CLAIM_REGISTER.csv`. Les statuts signifient :

- `PRESENT_VERIFIED` : le fichier ou la source existe et le contrôle indiqué a été réalisé.
- `PRESENT_UNVERIFIED` : l'affirmation est présente ou calculable partiellement, mais la preuve requise manque ou dépasse les sorties disponibles.
- `MISSING` : l'élément nécessaire n'a pas été trouvé.

Une affirmation numérique peut être `PRESENT_VERIFIED` comme calcul sur les fichiers fournis tout en restant insuffisante pour transformer une conjecture mathématique en théorème. C'est notamment le cas des 92 points.

Le registre doit être lu avec CORRECTIONS.md, LIMITS_AND_MISSING.md et les sorties de evidence/reanalysis/.

Dans la corrective v0.1a, C08 est volontairement borné au fait directement attesté : les deux constantes sont fixées pendant la réanalyse. Leur indépendance historique vis-à-vis du corpus n'est pas démontrée. C25 isole en outre la superellipse p=3 comme extension exploratoire hors du noyau à courbure strictement positive.
