#A faire

# Trouver les meilleures régressions  pour prédire hospital_los_day: length of stay in hospital (days, numeric)
# et  day_28_flg: death within 28 days (binary: 1 = yes, 0 = no)



data_tot<- readRDS("data.rds")

# informations sur la base :

### https://physionet.org/content/mimic2-iaccd/1.0/

## Chargement des données
library(ggplot2)
library(tidyverse)
library(class)




data_tot <- data_tot %>%
  select(hospital_los_day, day_28_flg, hosp_exp_flg, age, sofa_first, icu_los_day, 
         chf_flg, afib_flg, renal_flg, liver_flg, copd_flg, cad_flg, stroke_flg, 
         mal_flg, resp_flg, hr_1st, map_1st, temp_1st)

data_tot$hosp_exp_flg <- as.factor(data_tot$hosp_exp_flg)
set.seed(45)

sample_train <- sample(1:dim(data_tot)[1], 600, replace = FALSE)
data_train <- data_tot[sample_train, ]
data_test <- data_tot[-sample_train, ]

# Test de l'entraînement des modèles sur un fold
folds <- cut(seq(1, nrow(data_train)), breaks = 10, labels = FALSE)

# Modèle de régression logistique : prédiction de day_28_flg
fit <- glm(day_28_flg ~ hosp_exp_flg + age + sofa_first + icu_los_day + chf_flg + afib_flg + renal_flg + liver_flg + copd_flg + cad_flg + stroke_flg + mal_flg + resp_flg + hr_1st + map_1st + temp_1st, 
           data = data_train[!folds == 1, ], family = "binomial")
summary(fit)

predictions <- predict(fit, data_train[!folds == 1, ], type = "response")

# Seuil à 0.5
table(predictions > 0.5, data_train$day_28_flg[!folds == 1])

# Calcul des métriques de performance
conf_mat <- table(predictions > 0.5, data_train$day_28_flg[folds != 1])
TP <- conf_mat[2, 2]
FP <- conf_mat[2, 1]
FN <- conf_mat[1, 2]
TN <- conf_mat[1, 1]

precision <- TP / (TP + FP)
rappel <- TP / (TP + FN) # Sensibilité
specificite <- TN / (TN + FP)
vpn <- TN / (TN + FN)
acc <- (TP + TN) / sum(conf_mat)
F1 <- 2 * (precision * rappel) / (precision + rappel)

# Afficher les métriques
list(precision = precision, rappel = rappel, 
     specificite = specificite, vpn = vpn, acc = acc, F1 = F1)

# Cross-validation du meilleur seuil glm :
seuils <- seq(0.1, 0.9, by = 0.1)
res_data <- data.frame()
for(i in 1:10) { 
  for(seuil in seuils) {
    # Ajustement du modèle
    fit_i <- glm(day_28_flg ~ hosp_exp_flg + age + sofa_first + icu_los_day + chf_flg + afib_flg + renal_flg + liver_flg + copd_flg + cad_flg + stroke_flg + mal_flg + resp_flg + hr_1st + map_1st + temp_1st, 
                 data = data_train[folds != i, ], family = "binomial")
    predictions_i <- predict(fit_i, data_test, type = "response")
    
    # Matrice de confusion
    tableau_contingence <- table(predictions_i > seuil, data_test$day_28_flg)
    TP <- tableau_contingence[2, 2]
    FP <- tableau_contingence[2, 1]
    FN <- tableau_contingence[1, 2]
    TN <- tableau_contingence[1, 1]
    
    # Calcul des métriques
    precision <- TP / (TP + FP)
    rappel <- TP / (TP + FN)
    acc <- (TP + TN) / sum(tableau_contingence)
    F1 <- 2 * (precision * rappel) / (precision + rappel)
    
    # Stocker les résultats
    res_data <- rbind(res_data, 
                      data.frame(Fold = i, Seuil = seuil, Precision = precision, Rappel = rappel, 
                                 Accuracy = acc, F1 = F1))
  }
}

# Identifier le meilleur seuil
best_result <- res_data[which.max(res_data$F1), ]
print(best_result)

# Régression linéaire pour prédire hospital_los_day
set.seed(123)
sample_train <- sample(1:dim(data_tot)[1], 600, replace = FALSE)
data_train <- data_tot[sample_train, ]
data_test <- data_tot[-sample_train, ]

folds <- cut(seq(1, nrow(data_train)), breaks = 10, labels = FALSE)
res_data <- data.frame(Fold = integer(), MSE = numeric(), RMSE = numeric(), R2 = numeric())

for(i in 1:10) {
  train_fold <- data_train[folds != i, ]
  val_fold <- data_train[folds == i, ]
  
  fit <- lm(hospital_los_day ~ age + hosp_exp_flg + sofa_first + icu_los_day + chf_flg + afib_flg + renal_flg + liver_flg + copd_flg + cad_flg + stroke_flg + mal_flg + resp_flg + hr_1st + map_1st + temp_1st, 
            data = train_fold)
  
  predictions <- predict(fit, val_fold)
  
  actual <- val_fold$hospital_los_day
  mse <- mean((predictions - actual)^2)
  rmse <- sqrt(mse)
  r2 <- 1 - sum((predictions - actual)^2) / sum((actual - mean(actual))^2)
  
  res_data <- rbind(res_data, data.frame(Fold = i, MSE = mse, RMSE = rmse, R2 = r2))
}

mean_results <- colMeans(res_data[, -1])
print(mean_results)

final_model <- lm(hospital_los_day ~ age + hosp_exp_flg + sofa_first + icu_los_day + chf_flg + afib_flg + renal_flg + liver_flg + copd_flg + cad_flg + stroke_flg + mal_flg + resp_flg + hr_1st + map_1st + temp_1st, 
                  data = data_train)

final_predictions <- predict(final_model, data_test)
test_actual <- data_test$hospital_los_day

test_mse <- mean((final_predictions - test_actual)^2)
test_rmse <- sqrt(test_mse)
test_r2 <- 1 - sum((final_predictions - test_actual)^2) / sum((test_actual - mean(test_actual))^2)

list(Test_MSE = test_mse, Test_RMSE = test_rmse, Test_R2 = test_r2)

################################################################################################################
#################################################################################################
#Standarisation:
# Chargement des bibliothèques
library(ggplot2)
library(tidyverse)

# Charger les données
data_tot <- readRDS("data.rds")

# Sélectionner les colonnes pertinentes
data_tot <- data_tot %>% select(hospital_los_day, day_28_flg, hosp_exp_flg, age, sofa_first, icu_los_day, 
                                chf_flg, afib_flg, renal_flg, liver_flg, copd_flg, cad_flg, stroke_flg, 
                                mal_flg, resp_flg, hr_1st, map_1st, temp_1st)

# Standardiser les variables quantitatives
quantitative_vars <- c("age", "sofa_first", "icu_los_day", "hr_1st", "map_1st", "temp_1st")
data_tot[quantitative_vars] <- scale(data_tot[quantitative_vars])

# Transformer les variables catégorielles en facteurs
categorical_vars <- c("hosp_exp_flg", "chf_flg", "afib_flg", "renal_flg", "liver_flg", "copd_flg", 
                      "cad_flg", "stroke_flg", "mal_flg", "resp_flg")
data_tot[categorical_vars] <- lapply(data_tot[categorical_vars], as.factor)

# Diviser les données en train et test
set.seed(45)
sample_train <- sample(1:nrow(data_tot), 600, replace = FALSE)
data_train <- data_tot[sample_train, ]
data_test <- data_tot[-sample_train, ]

# Modèle de régression logistique pour prédire day_28_flg
fit_logistic <- glm(day_28_flg ~ ., data = data_train, family = "binomial")
summary(fit_logistic)

# Prédictions sur le jeu de test
predictions_logistic <- predict(fit_logistic, data_test, type = "response")

# Matrice de confusion pour évaluer le modèle
table(predictions_logistic > 0.5, data_test$day_28_flg)

# Modèle de régression linéaire pour prédire hospital_los_day
fit_linear <- lm(hospital_los_day ~ ., data = data_train)
summary(fit_linear)

# Prédictions sur le jeu de test
predictions_linear <- predict(fit_linear, data_test)

# Calcul des métriques pour le modèle linéaire
mse <- mean((predictions_linear - data_test$hospital_los_day)^2)
rmse <- sqrt(mse)
r2 <- 1 - sum((predictions_linear - data_test$hospital_los_day)^2) / sum((data_test$hospital_los_day - mean(data_test$hospital_los_day))^2)

# Afficher les résultats du modèle linéaire
list(MSE = mse, RMSE = rmse, R2 = r2)



