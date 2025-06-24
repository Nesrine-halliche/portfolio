###J'ai d'abond commencé à créer les modeles manuellement puis j'ai fait un autre code avec une liste et visualisation graphique:
#Charger les Library
library(tidymodels)
library(modeldata)

# Données
data(cells)
cells <- cells %>% select(-case)  
glimpse(cells)

# Séparation
set.seed(123)
cells_split <- initial_split(cells, strata = class)
cells_train <- training(cells_split)
cells_test <- testing(cells_split)

# Prétraitement
rec <- recipe(class ~ ., data = cells_train) %>%
  step_normalize(all_numeric_predictors())

# Modèle 1 : Régression Logistique
log_model <- logistic_reg() %>% set_engine("glm")
log_wf <- workflow() %>% add_model(log_model) %>% add_recipe(rec)
log_fit <- fit(log_wf, data = cells_train)
log_pred <- predict(log_fit, new_data = cells_test, type = "class")
log_metrics <- log_pred %>% bind_cols(cells_test) %>% metrics(truth = class, estimate = .pred_class)

# Modèle 2 : Arbre de décision
dt_model <- decision_tree() %>% set_engine("rpart") %>% set_mode("classification")
dt_wf <- workflow() %>% add_model(dt_model) %>% add_recipe(rec)
dt_fit <- fit(dt_wf, data = cells_train)
dt_pred <- predict(dt_fit, new_data = cells_test, type = "class")
dt_metrics <- dt_pred %>% bind_cols(cells_test) %>% metrics(truth = class, estimate = .pred_class)

# Modèle 3 : Random Forest
rf_model <- rand_forest() %>% set_engine("ranger") %>% set_mode("classification")
rf_wf <- workflow() %>% add_model(rf_model) %>% add_recipe(rec)
rf_fit <- fit(rf_wf, data = cells_train)
rf_pred <- predict(rf_fit, new_data = cells_test, type = "class")
rf_metrics <- rf_pred %>% bind_cols(cells_test) %>% metrics(truth = class, estimate = .pred_class)

# Modèle 4 : XGBoost
xgb_model <- boost_tree() %>% set_engine("xgboost") %>% set_mode("classification")
xgb_wf <- workflow() %>% add_model(xgb_model) %>% add_recipe(rec)
xgb_fit <- fit(xgb_wf, data = cells_train)
xgb_pred <- predict(xgb_fit, new_data = cells_test, type = "class")
xgb_metrics <- xgb_pred %>% bind_cols(cells_test) %>% metrics(truth = class, estimate = .pred_class)

# Comparer les performances
list(Logistic = log_metrics, DecisionTree = dt_metrics, RandomForest = rf_metrics, XGBoost = xgb_metrics)


###########################################################
##J'ai voulu ajouter une visualisation graphique donc j'ai fait un autre code:
library(ggplot2)
library(patchwork)

# Modèles
log_model <- logistic_reg() %>% set_engine("glm")
dt_model <- decision_tree() %>% set_engine("rpart") %>% set_mode("classification")
rf_model <- rand_forest() %>% set_engine("ranger") %>% set_mode("classification")
xgb_model <- boost_tree() %>% set_engine("xgboost") %>% set_mode("classification")

# Evaluer les modèles avec la fonction evaluate_model
evaluate_model <- function(model, name) {
  wf <- workflow() %>% add_model(model) %>% add_recipe(rec)
  fit <- fit(wf, data = cells_train)
  pred <- predict(fit, new_data = cells_test) %>%
    bind_cols(predict(fit, new_data = cells_test, type = "prob")) %>%
    bind_cols(cells_test) %>%
    mutate(Model = name)
  return(pred)
}

# Prédictions :
log_pred <- evaluate_model(log_model, "Logistic Regression")
dt_pred <- evaluate_model(dt_model, "Decision Tree")
rf_pred <- evaluate_model(rf_model, "Random Forest")
xgb_pred <- evaluate_model(xgb_model, "XGBoost")

# Combinaison
all_preds <- bind_rows(log_pred, dt_pred, rf_pred, xgb_pred)

# Matrice de Confusion
conf_matrix <- all_preds %>%
  conf_mat(truth = class, estimate = .pred_class) %>%
  autoplot(type = "heatmap") +  # Utilisation d'un type heatmap pour la matrice
  scale_fill_gradient(low = "white", high = "blue") +  # Gradient de couleurs
  labs(title = "Matrice de Confusion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Courbe ROC
roc_data <- all_preds %>%
  group_by(Model) %>%
  roc_curve(truth = class, .pred_PS)

roc_curve_plot <- ggplot(roc_data, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_line(size = 1.5) +  # Augmenter la taille de la ligne
  geom_abline(linetype = "dashed", color = "gray") +  # Ligne de référence pour l'égalité
  theme_minimal() +
  scale_color_brewer(palette = "Set1") +  # Palette de couleurs différente
  labs(title = "Courbes ROC", x = "1 - Spécificité", y = "Sensibilité")



metrics_plot <- all_preds %>%
  group_by(Model) %>%
  summarise(
    Accuracy = accuracy_vec(truth = class, estimate = .pred_class),
    AUC = roc_auc_vec(truth = class, estimate = .pred_PS)
  ) %>%
  pivot_longer(cols = c(Accuracy, AUC), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Value, y = Model, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +  # Barres horizontales
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") +  # Palette de couleurs différente
  labs(title = "Comparaison des performances des modèles", x = "Valeur", y = "Modèle")

#résulat
conf_matrix / roc_curve_plot / metrics_plot
