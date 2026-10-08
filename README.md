# Modélisation spatiale des expositions environnementales aux PFAS
### Étude de simulation — cas de la France

Ce projet R accompagne le rapport *« Modélisation spatiale des expositions
environnementales aux PFAS : revue de la littérature et étude de simulation
appliquée au cas de la France »*. Il reproduit intégralement l'étude de
simulation décrite dans la section 4 (méthodologie) et la section 5
(résultats) du rapport.

## Structure du projet

```
pfas_project/
├── R/
│   ├── 00_setup.R                    # packages, graine, parametres, hotspots
│   ├── 01_simulate_true_field.R      # simulation du champ "vrai" de concentration
│   ├── 02_sample_monitoring_network.R# echantillonnage du reseau de stations
│   ├── 03_variogram_kriging.R        # variogramme + krigeage ordinaire
│   ├── 04_cross_validation.R         # validation croisee (LOOCV), krigeage vs IDW
│   ├── 05_population_exposure.R      # exposition population + simulation conditionnelle
│   └── 06_run_all.R                  # script maitre (execute tout, dans l'ordre)
├── data/                             # objets .rds generes (champ vrai, stations, predictions...)
└── output/
    ├── figures/                      # cartes et graphiques (.png)
    └── tables/                       # tableaux de resultats (.csv)
```

## Prérequis

R (>= 4.1) et les packages suivants :

```r
install.packages(c("sf", "sp", "gstat", "ggplot2", "viridis", "gridExtra"))
```

Sous Ubuntu/Debian, ces packages peuvent aussi être installés via apt :

```bash
sudo apt-get install r-cran-sf r-cran-sp r-cran-gstat r-cran-ggplot2 \
                      r-cran-viridis r-cran-gridextra
```

## Reproduire l'étude

Depuis la racine du projet (`pfas_project/`) :

```bash
Rscript R/06_run_all.R
```

Ce script exécute successivement les six scripts numérotés et régénère
l'ensemble des fichiers de `data/`, `output/figures/` et `output/tables/`.
La graine aléatoire étant fixée (`set.seed(20260824)` dans `00_setup.R`),
les résultats sont strictement reproductibles.

Chaque script peut aussi être exécuté séparément (dans l'ordre), à
condition d'être lancé depuis la racine du projet, par exemple :

```bash
Rscript R/01_simulate_true_field.R
Rscript R/02_sample_monitoring_network.R
Rscript R/03_variogram_kriging.R
Rscript R/04_cross_validation.R
Rscript R/05_population_exposure.R
```

## Résultats principaux attendus

- **Validation croisée** (`output/tables/table_cv_metrics.csv`) : le
  krigeage ordinaire (RMSE ≈ 8,1 ng/L, R² ≈ 0,66) surpasse l'interpolation
  par distance inverse — IDW — (RMSE ≈ 10,4 ng/L, R² ≈ 0,52).
- **Exposition de la population**
  (`output/tables/table_exposition_population.csv`) : l'estimation
  ponctuelle de la population exposée à un dépassement du seuil
  réglementaire européen (100 ng/L) à partir de la seule carte krigée
  sous-estime fortement la population réellement exposée dans le champ
  simulé ; les simulations géostatistiques conditionnelles permettent de
  quantifier cette incertitude sous forme d'un intervalle de crédibilité.

## Adapter le projet à des données réelles

Pour appliquer cette chaîne méthodologique aux données réelles publiées
par le BRGM (plateforme PFAS France, 2025, > 2,3 millions d'analyses), il
suffit de remplacer, dans `02_sample_monitoring_network.R`, l'objet
`stations` simulé par un import des données réelles (coordonnées,
concentration mesurée, indicateur de censure), puis de relancer les
scripts `03` à `05` sans autre modification. La grille de calcul et la
grille de population (`00_setup.R` et `05_population_exposure.R`)
devraient alors être remplacées par des données géographiques et
démographiques officielles (IGN, INSEE, Eurostat) pour une étude
appliquée réelle plutôt qu'une simulation pédagogique.
