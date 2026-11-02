# 🌡️ Détection de ruptures dans la série de températures CET

Analyse de la série **CET (Central England Temperature)** du Met Office —
la plus longue série de températures instrumentales au monde, de 1659 à
aujourd'hui — pour y détecter des ruptures : des moments où la température
moyenne, ou la tendance du réchauffement, change brutalement.

Deux approches sont comparées :
- une **programmation dynamique exacte**, qui minimise un coût pénalisé (BIC)
- un **algorithme génétique**, plus rapide pour explorer de grands espaces de
  ruptures, avec une variante tenant compte d'erreurs autocorrélées (AR(1))

## Aperçu

- 📈 Visualisation de la série brute, moyenne mobile et bande de confiance
- 📉 Régression linéaire globale et série détendancée
- 🔍 Détection de ruptures **en moyenne** (changement de niveau)
- 🔍 Détection de ruptures **en tendance** (changement de pente)
- 🧬 Version algo génétique des deux détections, y compris avec erreurs AR(1)
- ✅ Vérification de l'autocorrélation des résidus

## Données

`meanstemp.xlsx` — températures moyennes mensuelles et annuelles (°C),
1659–2021, telles que publiées par le Met Office.

| Colonne | Contenu |
|---|---|
| `Annee` | année |
| `JAN`…`DEC` | température moyenne du mois |
| `YEAR` | température moyenne annuelle |

La dernière ligne (2021, encore incomplète au moment de la collecte) est
retirée automatiquement par le script.

## Installation

```r
install.packages(c(
  "readxl", "ggplot2", "forecast", "ggfortify", "tidyverse",
  "timeSeries", "xts", "zoo", "GA", "microbenchmark"
))
```

## Utilisation

1. Placez `meanstemp.xlsx` et `memoire.R` dans le même dossier.
2. Ouvrez `memoire.R` dans RStudio et exécutez-le section par section — le
   script est pensé comme un carnet d'exploration, pas comme un pipeline à
   lancer d'un bloc.

Grandes sections, dans l'ordre :

```
Chargement & nettoyage des données
        │
        ├── Visualisation (série brute, moyenne mobile, bande de confiance)
        ├── Régression linéaire globale + série détendancée
        │
        ├── Ruptures en MOYENNE  ──┐
        ├── Ruptures en TENDANCE ──┼── par programmation dynamique + BIC
        │                          │
        ├── Autocorrélation des résidus
        │
        ├── Ruptures en MOYENNE  ──┐
        ├── Ruptures en TENDANCE ──┼── par algorithme génétique
        └── Variante avec erreurs AR(1) ──┘
```

## ⚠️ Bon à savoir

- Le calcul du coût pour les ruptures **en tendance** (`cout_t`) refait une
  régression linéaire pour chaque paire de points de la série : comptez
  quelques minutes sur la série complète (362 points).
- Avec un grand nombre de ruptures autorisées (`Kmax`), le critère BIC a
  tendance à **sur-détecter** — n'hésitez pas à regarder la courbe des coûts
  pénalisés en fonction de k pour repérer le "coude" plutôt que de prendre
  le minimum brut.
- Les runs de l'algorithme génétique (`popSize` 250-500, `maxiter` 300-400)
  peuvent aussi prendre plusieurs minutes chacun.

## Historique des corrections

- Chemin de fichier en dur remplacé par un chemin relatif, pour que le script
  tourne sur n'importe quelle machine.
- Bug dans la droite de régression tracée sur la série (`f()`) : intercept et
  pente étaient inversés, ce qui faisait diverger la droite de manière
  absurde en fin de série.
- Bug dans la sélection du nombre optimal de ruptures en tendance
  (`biccorr()`) : mauvaise syntaxe d'indexation de liste, corrigée.
- Bug de portée dans une fonction auxiliaire de l'algorithme génétique
  (`rec()`), qui référençait une variable inexistante.
- Nettoyage : suppression d'une dépendance chargée mais jamais utilisée.

## Licence

Libre d'utilisation pour un usage académique ou personnel.
