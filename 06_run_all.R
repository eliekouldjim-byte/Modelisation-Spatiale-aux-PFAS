## ============================================================
## 06_run_all.R
## Script maitre : execute l'ensemble du pipeline dans l'ordre.
## Usage : Rscript R/06_run_all.R   (a lancer depuis la racine
## du projet, pfas_project/)
## ============================================================

scripts <- c(
  "R/00_setup.R",
  "R/01_simulate_true_field.R",
  "R/02_sample_monitoring_network.R",
  "R/03_variogram_kriging.R",
  "R/04_cross_validation.R",
  "R/05_population_exposure.R"
)

for (s in scripts) {
  cat("\n============================================================\n")
  cat("Execution de :", s, "\n")
  cat("============================================================\n")
  source(s)
}

cat("\nPipeline complet termine. Resultats dans data/, output/figures/, output/tables/.\n")
