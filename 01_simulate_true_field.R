## ============================================================
## 01_simulate_true_field.R
## Simulation du champ "vrai" de concentration en PFAS
## (somme PFAS, en ng/L) sur une grille reguliere couvrant le pays.
## ============================================================
## Logique de simulation :
##   1) Tendance deterministe = decroissance exponentielle de la
##      concentration autour de "hotspots" industriels/urbains,
##      plus un fond regional bas (bruit de fond diffus, coherent
##      avec les niveaux nationaux rapportes par Anses 2023-2025).
##   2) Processus spatial aleatoire correle (champ gaussien) simule
##      par methode de bandes tournantes via gstat::gstatSim /
##      krige non conditionnel, pour representer l'heterogeneite
##      non expliquee par la tendance (variabilite locale des sols,
##      hydrogeologie, rejets diffus).
##   3) Le champ simule sert de "verite terrain" (population de
##      reference) a partir de laquelle on tire un echantillon de
##      stations de mesure (etape 02), afin d'evaluer la capacite
##      des modeles geostatistiques a la reconstruire (etape 03-04).
## ============================================================

source("R/00_setup.R")

## ---- 1. Grille reguliere sur l'emprise du pays -------------------------
n_side <- 80   # grille n_side x n_side (resolution ~ 12 km)
grid_xy <- expand.grid(
  lon = seq(bbox_France["xmin"], bbox_France["xmax"], length.out = n_side),
  lat = seq(bbox_France["ymin"], bbox_France["ymax"], length.out = n_side)
)

## Masque grossier "en forme de France" : on retire les coins trop
## eloignes de tout hotspot et de la diagonale principale du territoire,
## pour eviter un pave parfaitement rectangulaire (simplification
## acceptable pour une etude de simulation pedagogique).
dist_to_center <- sqrt((grid_xy$lon - 2.5)^2 / 3.2^2 + (grid_xy$lat - 46.5)^2 / 3.5^2)
grid_xy <- grid_xy[dist_to_center < 1.35, ]
rownames(grid_xy) <- NULL

## ---- 2. Tendance deterministe : sources + fond regional ----------------
## Distance (en degres, approx. 1 deg ~ 100 km) a chaque hotspot,
## contribution en decroissance exponentielle (portee ~ 60 km).
portee_km   <- 35
deg_per_km  <- 1 / 100

trend <- rep(0, nrow(grid_xy))
for (i in seq_len(nrow(hotspots))) {
  d_km <- sqrt((grid_xy$lon - hotspots$lon[i])^2 +
               (grid_xy$lat - hotspots$lat[i])^2) / deg_per_km
  trend <- trend + hotspots$intensity[i] * exp(-d_km / portee_km)
}

fond_regional <- 3.0    # ng/L, bruit de fond diffus (TFA/PFAS ultra-courts, cf. Anses)
concentration_max_source <- 130  # ng/L, contribution max pres d'une source forte

trend_ngL <- fond_regional + concentration_max_source * trend

## ---- 3. Composante spatiale aleatoire correlee (bruit structure) -------
coordinates(grid_xy) <- ~lon + lat
proj4string(grid_xy) <- CRS("+proj=longlat +datum=WGS84")

## Modele de variogram "non conditionnel" utilise pour simuler un champ
## gaussien correle (moyenne nulle, portee courte = heterogeneite locale).
vgm_sim <- vgm(psill = 9, model = "Exp", range = 0.45, nugget = 1.5)
g_sim <- gstat(formula = z ~ 1, locations = grid_xy, dummy = TRUE,
               beta = 0, model = vgm_sim, nmax = 40)
sim <- predict(g_sim, newdata = grid_xy, nsim = 1, debug.level = 0)

## ---- 4. Champ vrai final : tendance + bruit spatial + plancher a 0 ------
true_field <- data.frame(
  lon = coordinates(grid_xy)[, 1],
  lat = coordinates(grid_xy)[, 2],
  trend_ngL = trend_ngL,
  noise_ngL = sim$sim1,
  true_ngL  = pmax(0, trend_ngL + sim$sim1)
)

saveRDS(true_field, file.path(dir_data, "true_field.rds"))

cat(sprintf("Champ vrai simule sur %d points de grille.\n", nrow(true_field)))
cat(sprintf("  - concentration moyenne  : %.2f ng/L\n", mean(true_field$true_ngL)))
cat(sprintf("  - concentration mediane  : %.2f ng/L\n", median(true_field$true_ngL)))
cat(sprintf("  - concentration max      : %.2f ng/L\n", max(true_field$true_ngL)))
cat(sprintf("  - %% de points > 100 ng/L (seuil UE somme PFAS) : %.1f%%\n",
            100 * mean(true_field$true_ngL > 100)))
