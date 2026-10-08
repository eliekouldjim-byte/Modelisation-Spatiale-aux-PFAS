## ============================================================
## 00_setup.R
## Modelisation spatiale des expositions environnementales aux PFAS
## Etude de simulation - cas de la France
## ============================================================
## Charge les packages necessaires et definit les parametres
## globaux (graine aleatoire, dossiers de sortie, emprise
## geographique du "pays" simule).
## ============================================================

pkgs <- c("sf", "sp", "gstat", "ggplot2", "viridis")
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(20260824)   # graine pour reproductibilite

## Dossiers
dir_root   <- "/home/claude/pfas_project"
dir_data   <- file.path(dir_root, "data")
dir_fig    <- file.path(dir_root, "output/figures")
dir_tab    <- file.path(dir_root, "output/tables")
dir.create(dir_data, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_fig,  showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tab,  showWarnings = FALSE, recursive = TRUE)

## Emprise approximative de la France metropolitaine (boite englobante,
## en degres decimaux, CRS = EPSG:4326)
bbox_France <- c(xmin = -4.9, xmax = 8.3, ymin = 41.2, ymax = 51.2)

## Points "sources industrielles / hotspots" utilises pour generer
## un champ de concentration PFAS realiste (inspires de la litterature :
## agglomeration lyonnaise, region parisienne, facade Nord (Dunkerque),
## Sud-Ouest -- cf. BRGM 2025 ; Munoz et al. 2025).
hotspots <- data.frame(
  name = c("Lyon", "Paris", "Dunkerque", "Sud-Ouest", "Marseille"),
  lon  = c(4.83,   2.35,    2.38,        -0.5,        5.37),
  lat  = c(45.76,  48.85,   51.03,       44.5,        43.30),
  intensity = c(1.60, 0.55, 0.70, 0.45, 0.50)   # poids relatif de la source
)

cat("Setup termine : packages charges, dossiers crees, parametres definis.\n")
