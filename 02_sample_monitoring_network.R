## ============================================================
## 02_sample_monitoring_network.R
## Simulation d'un reseau d'observation (stations de prelevement)
## a partir du champ vrai, avec:
##   - un echantillonnage preferentiel (plus dense pres des zones
##     urbaines/industrielles, comme le reseau reel BRGM/Anses,
##     cf. BRGM 2025 : 21 000 sites, forte densite pres des
##     agglomerations et sites industriels) ;
##   - une erreur de mesure additive (incertitude analytique) ;
##   - une censure a gauche sous la limite de quantification (LOQ),
##     frequente pour les composes traces (cf. Smalling et al. 2023).
## ============================================================

source("R/00_setup.R")
true_field <- readRDS(file.path(dir_data, "true_field.rds"))

n_stations <- 220   # ordre de grandeur d'un reseau de surveillance regional

## ---- Probabilite d'echantillonnage preferentiel -------------------------
## Plus forte pres des hotspots (comme dans la realite : la surveillance
## est renforcee autour des sites suspects) + une composante uniforme pour
## représenter le reseau de fond (ex-ante, non informe par la contamination).
w_pref <- 0
for (i in seq_len(nrow(hotspots))) {
  d_km <- sqrt((true_field$lon - hotspots$lon[i])^2 +
               (true_field$lat - hotspots$lat[i])^2) * 100
  w_pref <- w_pref + exp(-d_km / 80)
}
w_unif <- 1
w_total <- 0.55 * w_pref / max(w_pref) + 0.45 * w_unif

set.seed(2026)
idx_sample <- sample(seq_len(nrow(true_field)), size = n_stations,
                      prob = w_total, replace = FALSE)

stations <- true_field[idx_sample, c("lon", "lat", "true_ngL")]

## Leger "jitter" des coordonnees (+/- ~2 km) : les stations reelles ne sont
## jamais exactement au centre des mailles d'une grille de calcul ; ce
## decalage evite aussi toute coincidence exacte entre points de mesure et
## noeuds de la grille de prediction (necessaire pour la simulation
## geostatistique conditionnelle a l'etape 05).
stations$lon <- stations$lon + rnorm(n_stations, 0, 0.01)
stations$lat <- stations$lat + rnorm(n_stations, 0, 0.01)

## ---- Erreur de mesure et censure sous la LOQ -----------------------------
sigma_mesure <- 1.8   # ng/L, ecart-type de l'erreur analytique
LOQ <- 2               # ng/L, limite de quantification usuelle pour la somme PFAS

stations$mesure_ngL <- pmax(0, stations$true_ngL + rnorm(n_stations, 0, sigma_mesure))
stations$censure    <- stations$mesure_ngL < LOQ
stations$mesure_ngL_analyse <- ifelse(stations$censure, LOQ / 2, stations$mesure_ngL)

saveRDS(stations, file.path(dir_data, "stations_echantillon.rds"))
write.csv(stations, file.path(dir_tab, "stations_echantillon.csv"), row.names = FALSE)

cat(sprintf("Reseau simule : %d stations.\n", nrow(stations)))
cat(sprintf("  - %% de mesures censurees (< LOQ = %.0f ng/L) : %.1f%%\n",
            LOQ, 100 * mean(stations$censure)))
cat(sprintf("  - concentration mesuree moyenne : %.2f ng/L\n", mean(stations$mesure_ngL_analyse)))
