## ============================================================
## 03_variogram_kriging.R
## Estimation du variogramme empirique, ajustement d'un modele
## theorique, et krigeage ordinaire pour reconstruire la surface
## de concentration en PFAS a partir du reseau de stations simule.
## ============================================================

source("R/00_setup.R")
true_field <- readRDS(file.path(dir_data, "true_field.rds"))
stations   <- readRDS(file.path(dir_data, "stations_echantillon.rds"))

## ---- Objets spatiaux (sp) ------------------------------------------------
stations_sp <- stations
coordinates(stations_sp) <- ~lon + lat
proj4string(stations_sp) <- CRS("+proj=longlat +datum=WGS84")

grid_sp <- true_field[, c("lon", "lat")]
coordinates(grid_sp) <- ~lon + lat
proj4string(grid_sp) <- CRS("+proj=longlat +datum=WGS84")
gridded(grid_sp) <- FALSE  # points irreguliers -> pas de grille sp stricte

## ---- Variogramme empirique et ajustement ---------------------------------
## Remarque : les coordonnees sont en longitude/latitude (CRS geographique) ;
## gstat calcule alors automatiquement les distances en kilometres
## (grand cercle), d'ou des parametres cutoff/width exprimes en km.
## On travaille sur la concentration mesuree (LOQ/2 pour les valeurs censurees).
v_emp <- variogram(mesure_ngL_analyse ~ 1, data = stations_sp, cutoff = 400, width = 25)

## Valeurs initiales raisonnables (nugget/sill/portee, en km) puis ajustement
## automatique par moindres carres ponderes (fit.variogram).
v_ini <- vgm(psill = var(stations$mesure_ngL_analyse), model = "Exp",
             range = 60, nugget = 3)
v_fit <- fit.variogram(v_emp, model = v_ini)

cat("Modele de variogramme ajuste :\n")
print(v_fit)

## ---- Graphique du variogramme --------------------------------------------
png(file.path(dir_fig, "fig_variogram.png"), width = 1400, height = 1000, res = 180)
plot(v_emp, v_fit, main = "Variogramme empirique et modele ajuste (somme PFAS)",
     xlab = "Distance (km)", ylab = "Semi-variance")
dev.off()

## ---- Krigeage ordinaire sur la grille -------------------------------------
ok_pred <- krige(mesure_ngL_analyse ~ 1, locations = stations_sp,
                  newdata = grid_sp, model = v_fit)

pred_df <- data.frame(
  lon = coordinates(grid_sp)[, 1],
  lat = coordinates(grid_sp)[, 2],
  pred_ngL = ok_pred$var1.pred,
  var_pred = ok_pred$var1.var
)
pred_df$pred_ngL <- pmax(0, pred_df$pred_ngL)
pred_df$true_ngL <- true_field$true_ngL

saveRDS(pred_df, file.path(dir_data, "kriging_predictions.rds"))
saveRDS(v_fit,   file.path(dir_data, "variogram_model.rds"))

cat(sprintf("\nKrigeage effectue sur %d points de grille.\n", nrow(pred_df)))
cat(sprintf("  - concentration predite moyenne : %.2f ng/L (vrai : %.2f ng/L)\n",
            mean(pred_df$pred_ngL), mean(pred_df$true_ngL)))

## ---- Cartes : champ vrai vs. champ krige ----------------------------------
p1 <- ggplot(true_field, aes(lon, lat, fill = true_ngL)) +
  geom_tile() +
  scale_fill_viridis(name = "ng/L", option = "inferno", limits = c(0, 100), oob = scales::squish) +
  coord_equal() + theme_minimal() +
  labs(title = "Champ 'vrai' simule - somme PFAS", x = "Longitude", y = "Latitude")

p2 <- ggplot(pred_df, aes(lon, lat, fill = pred_ngL)) +
  geom_tile() +
  geom_point(data = stations, aes(lon, lat), inherit.aes = FALSE,
             shape = 21, size = 0.6, color = "white", stroke = 0.2, alpha = 0.6) +
  scale_fill_viridis(name = "ng/L", option = "inferno", limits = c(0, 100), oob = scales::squish) +
  coord_equal() + theme_minimal() +
  labs(title = "Surface krigee (krigeage ordinaire)", x = "Longitude", y = "Latitude",
       subtitle = "Points blancs = stations du reseau simule (n = 220)")

ggsave(file.path(dir_fig, "fig_true_vs_kriged.png"),
       gridExtra::arrangeGrob(p1, p2, ncol = 2), width = 12, height = 5.2, dpi = 200)

p3 <- ggplot(pred_df, aes(lon, lat, fill = sqrt(var_pred))) +
  geom_tile() +
  scale_fill_viridis(name = "Ecart-type\nde krigeage\n(ng/L)", option = "viridis") +
  coord_equal() + theme_minimal() +
  labs(title = "Incertitude de la prediction (ecart-type de krigeage)",
       x = "Longitude", y = "Latitude")
ggsave(file.path(dir_fig, "fig_kriging_uncertainty.png"), p3, width = 6.5, height = 5.2, dpi = 200)

cat("Figures enregistrees dans output/figures/\n")
