## ============================================================
## 05_population_exposure.R
## Couplage de la surface krigee avec une grille de population
## simulee (densite plus forte pres des grandes agglomerations)
## pour estimer :
##   - l'exposition moyenne ponderee par la population ;
##   - la population estimee exposee au-dessus du seuil
##     reglementaire europeen (100 ng/L, somme de 20 PFAS,
##     directive UE 2020/2184) ;
##   - une carte de risque de depassement (P[concentration > seuil]),
##     obtenue par simulation conditionnelle geostatistique afin de
##     propager l'incertitude du krigeage (cf. Goovaerts, 1997).
## ============================================================

source("R/00_setup.R")
true_field <- readRDS(file.path(dir_data, "true_field.rds"))
pred_df    <- readRDS(file.path(dir_data, "kriging_predictions.rds"))
stations   <- readRDS(file.path(dir_data, "stations_echantillon.rds"))
v_fit      <- readRDS(file.path(dir_data, "variogram_model.rds"))

SEUIL_UE <- 100  # ng/L, seuil reglementaire "somme PFAS" (directive UE 2020/2184)

## ---- 1. Grille de population simulee --------------------------------------
## Densite de population decroissant avec la distance aux grands poles
## urbains (memes centres que les hotspots "Paris" et "Lyon" + Marseille),
## calibree pour donner un ordre de grandeur national plausible.
pop_centers <- hotspots[hotspots$name %in% c("Paris", "Lyon", "Marseille"), ]
pop_centers$poids_pop <- c(9, 3, 2.5)  # poids demographique relatif approximatif

dens <- rep(0.15, nrow(true_field))  # densite de fond (zones rurales)
for (i in seq_len(nrow(pop_centers))) {
  d_km <- sqrt((true_field$lon - pop_centers$lon[i])^2 +
               (true_field$lat - pop_centers$lat[i])^2) * 100
  dens <- dens + pop_centers$poids_pop[i] * exp(-d_km / 40)
}

## Population totale simulee calee approximativement sur la France (~68 millions)
pop_totale_cible <- 68e6
poids_norm <- dens / sum(dens)
true_field$population <- poids_norm * pop_totale_cible
pred_df$population <- true_field$population

## ---- 2. Exposition moyenne ponderee par la population ---------------------
exp_vraie_pond   <- weighted.mean(true_field$true_ngL, true_field$population)
exp_predite_pond <- weighted.mean(pred_df$pred_ngL,   pred_df$population)

pop_exposee_vraie   <- sum(true_field$population[true_field$true_ngL > SEUIL_UE])
pop_exposee_predite <- sum(pred_df$population[pred_df$pred_ngL > SEUIL_UE])

cat(sprintf("Exposition moyenne ponderee par la population :\n"))
cat(sprintf("  - champ vrai      : %.2f ng/L\n", exp_vraie_pond))
cat(sprintf("  - champ krige     : %.2f ng/L\n", exp_predite_pond))
cat(sprintf("\nPopulation estimee exposee a > %d ng/L :\n", SEUIL_UE))
cat(sprintf("  - d'apres le champ vrai  : %s habitants\n", format(round(pop_exposee_vraie), big.mark = " ")))
cat(sprintf("  - d'apres le champ krige : %s habitants\n", format(round(pop_exposee_predite), big.mark = " ")))

## ---- 3. Simulations geostatistiques conditionnelles (propagation d'incertitude) ----
stations_sp <- stations
coordinates(stations_sp) <- ~lon + lat
proj4string(stations_sp) <- CRS("+proj=longlat +datum=WGS84")

grid_sp <- true_field[, c("lon", "lat")]
coordinates(grid_sp) <- ~lon + lat
proj4string(grid_sp) <- CRS("+proj=longlat +datum=WGS84")

## Un nugget numeriquement nul peut rendre le systeme de krigeage simule
## instable (matrice de covariance quasi singuliere) ; on impose un
## nugget plancher minime (stabilisation numerique standard, sans
## consequence pratique sur les resultats).
v_fit_sim <- v_fit
v_fit_sim$psill[v_fit_sim$model == "Nug"] <- max(0.5, v_fit_sim$psill[v_fit_sim$model == "Nug"])

n_sim <- 100
cond_sim <- krige(mesure_ngL_analyse ~ 1, locations = stations_sp, newdata = grid_sp,
                   model = v_fit_sim, nsim = n_sim, nmax = 40, debug.level = 0)

sim_mat <- as.data.frame(cond_sim)
sim_cols <- grep("^sim", names(sim_mat), value = TRUE)
sim_mat[sim_cols] <- lapply(sim_mat[sim_cols], function(x) pmax(0, x))

## Probabilite (frequence empirique) de depassement du seuil en chaque point
proba_depassement <- rowMeans(sim_mat[sim_cols] > SEUIL_UE)
pred_df$proba_depassement <- proba_depassement

## Distribution simulee de la population totale exposee (incertitude)
pop_exp_sim <- sapply(sim_cols, function(cn) {
  sum(true_field$population[sim_mat[[cn]] > SEUIL_UE])
})
ic95 <- quantile(pop_exp_sim, c(0.025, 0.5, 0.975))
cat(sprintf("\nIntervalle de credibilite (simulations conditionnelles, n=%d) pour la\n", n_sim))
cat(sprintf("population exposee a > %d ng/L :\n", SEUIL_UE))
cat(sprintf("  - mediane : %s habitants\n", format(round(ic95[2]), big.mark = " ")))
cat(sprintf("  - IC95%%   : [%s ; %s] habitants\n",
            format(round(ic95[1]), big.mark = " "), format(round(ic95[3]), big.mark = " ")))

saveRDS(pred_df, file.path(dir_data, "kriging_predictions.rds"))
write.csv(data.frame(indicateur = c("Exposition moyenne ponderee - vrai",
                                     "Exposition moyenne ponderee - krige",
                                     "Population exposee (vrai)",
                                     "Population exposee (krige)",
                                     "Population exposee - mediane simulee",
                                     "Population exposee - borne IC95 basse",
                                     "Population exposee - borne IC95 haute"),
                      valeur = c(exp_vraie_pond, exp_predite_pond,
                                 pop_exposee_vraie, pop_exposee_predite,
                                 ic95[2], ic95[1], ic95[3])),
          file.path(dir_tab, "table_exposition_population.csv"), row.names = FALSE)

## ---- 4. Cartes finales -----------------------------------------------------
p_risk <- ggplot(pred_df, aes(lon, lat, fill = proba_depassement)) +
  geom_tile() +
  scale_fill_viridis(name = "P(somme PFAS\n> 100 ng/L)", option = "magma", limits = c(0, 1)) +
  coord_equal() + theme_minimal() +
  labs(title = "Probabilite de depassement du seuil reglementaire UE",
       subtitle = sprintf("Estimee par %d simulations geostatistiques conditionnelles", n_sim),
       x = "Longitude", y = "Latitude")
ggsave(file.path(dir_fig, "fig_risk_map.png"), p_risk, width = 7, height = 5.5, dpi = 200)

p_pop <- ggplot(true_field, aes(lon, lat, fill = population)) +
  geom_tile() +
  scale_fill_viridis(name = "Population\n(hab./maille)", option = "cividis", trans = "sqrt") +
  coord_equal() + theme_minimal() +
  labs(title = "Grille de population simulee", x = "Longitude", y = "Latitude")
ggsave(file.path(dir_fig, "fig_population_grid.png"), p_pop, width = 6.5, height = 5.2, dpi = 200)

cat("\nCartes de risque et de population enregistrees dans output/figures/\n")
