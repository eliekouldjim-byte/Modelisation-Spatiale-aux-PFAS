## ============================================================
## 04_cross_validation.R
## Validation croisee "leave-one-out" (LOOCV) du modele de
## krigeage ordinaire, et comparaison a une methode naive
## (interpolation par distance inverse, IDW) utilisee comme
## reference de base (benchmark).
## ============================================================

source("R/00_setup.R")
stations <- readRDS(file.path(dir_data, "stations_echantillon.rds"))
v_fit    <- readRDS(file.path(dir_data, "variogram_model.rds"))

stations_sp <- stations
coordinates(stations_sp) <- ~lon + lat
proj4string(stations_sp) <- CRS("+proj=longlat +datum=WGS84")

## ---- LOOCV krigeage ordinaire ---------------------------------------------
cv_ok <- krige.cv(mesure_ngL_analyse ~ 1, locations = stations_sp, model = v_fit, nfold = nrow(stations))

## ---- LOOCV IDW (puissance 2) - methode de reference -----------------------
cv_idw <- krige.cv(mesure_ngL_analyse ~ 1, locations = stations_sp, set = list(idp = 2), nfold = nrow(stations))

metrics <- function(obs, pred) {
  err <- obs - pred
  c(ME   = mean(err),
    MAE  = mean(abs(err)),
    RMSE = sqrt(mean(err^2)),
    R2   = cor(obs, pred)^2)
}

res_ok  <- metrics(cv_ok$observed,  cv_ok$observed - cv_ok$residual)
res_idw <- metrics(cv_idw$observed, cv_idw$observed - cv_idw$residual)

tab_cv <- rbind(Krigeage_ordinaire = res_ok, IDW_puissance2 = res_idw)
print(round(tab_cv, 3))
write.csv(round(tab_cv, 3), file.path(dir_tab, "table_cv_metrics.csv"))

## ---- Graphique observe vs. predit (krigeage) ------------------------------
df_cv <- data.frame(observed = cv_ok$observed, predicted = cv_ok$observed - cv_ok$residual)
p_cv <- ggplot(df_cv, aes(observed, predicted)) +
  geom_point(alpha = 0.6, color = "#3b528b") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Validation croisee (leave-one-out) - krigeage ordinaire",
       x = "Concentration observee (ng/L)", y = "Concentration predite (ng/L)")
ggsave(file.path(dir_fig, "fig_cv_scatter.png"), p_cv, width = 6, height = 5.2, dpi = 200)

cat("\nValidation croisee terminee. Resultats enregistres dans output/tables/table_cv_metrics.csv\n")
