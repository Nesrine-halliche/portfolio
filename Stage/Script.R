choose.files()
library(dplyr)
library(summarytools)
library(skimr)
library(WeightIt)
library(cobalt)
library(broom)
library(ggplot2)
# 1. Charger les données
df <- read.csv("C:\\Users\\PC PEO MAX\\Downloads\\CRFMarianeDataSynthetiques_TVAE.csv", sep = ",", header = TRUE)

# 2.Aperçu complet de toutes les variables
dfSummary(df)
# Créer le résumé
tableau <- dfSummary(df)
view(tableau)


# Générer un résumé statistique propre avec skimr
descriptif <- skim(df)
View(descriptif)

#Tableau pour les valeurs manquantes :NA
na_table <- data.frame(
  variable = names(df),
  nb_NA = sapply(df, function(x) sum(is.na(x)))
) %>% filter(nb_NA > 0)

View(na_table)


# Variables explicatives 
vars <- c(
  "age", "sexe", "diabete", "hta", "tabacActif", "obesiteSurpoids",
  "atcdCardioVasculaire1", "atcdRespiratoires", "scoreKnaus", "glasgowIoa",
  "sous02Ioa", "fc", "tAs", "tAd", "temperature", 
  "moarriveeUrgences", "delaiDepuisApparitionPremierSymptome",
  "consultPrealablePassageSu", "presentationCliniqueUrgences",
  "creatinineAdmission", "dfg", "natremie", "crpAdmission", "leucocytesTotadmission"
)
df_clean <- df %>%
  mutate(across(all_of(vars), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

# 5. Sélection des colonnes utiles
df_ipw <- df_clean %>%
  select(A, Y, all_of(vars))
               #####################################
######################################################################
#Analyse brute (sans ajustement)
# 1. Modèle logistique simple
modele_brut <- glm(Y ~ A, data = df_ipw, family = binomial)# Charger les bons packages
# Générer les résultats OR + IC95%
effet_brut <- tidy(modele_brut, exponentiate = TRUE, conf.int = TRUE)

# Filtrer la variable A et formater le tableau
effet_brut_A <- effet_brut %>%
  filter(term == "A") %>%
  mutate(
    Variable = "A (antibio)",
    OR = round(estimate, 2),
    `IC 95%` = paste0("[", round(conf.low, 2), " ; ", round(conf.high, 2), "]"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  select(Variable, OR, `IC 95%`, `p-value`)

print(effet_brut_A)

               ####################################
###################################################################
#Estimation des poids IPTW (par score de propension)
w.out <- weightit(A ~ ., data = df_ipw, method = "ps", estimand = "ATE")

#Vérification du déséquilibre avant/après pondération
tab_unw <- bal.tab(w.out, un = TRUE, m.threshold = 0.1, disp.v.ratio = TRUE)

# Graphique de déséquilibre
love.plot(w.out, var.order = "unadjusted", threshold = 0.1,
          abs = TRUE, colors = c("red", "blue"), title = "Déséquilibre avant/après pondération")

#Affichage du tableau de déséquilibre
print(tab_unw)

#Ajout des poids au dataset
df_ipw$poids <- get.w(w.out)

#Régression pondérée pour effet marginal du traitement
df_ipw$Y <- factor(df_ipw$Y, levels = c(0, 1))  # Important pour binomial

modele_iptw <- glm(Y ~ A, data = df_ipw, weights = poids, family = binomial)

#Résumé du modèle
summary(modele_iptw)

#Coefficients exponentiés (OR) avec IC
tidy(modele_iptw, exponentiate = TRUE, conf.int = TRUE)

############################################################################
####################################entropy balancing##########################

#Étape 1 : Sélection des covariables pertinentes 
vars_ebal <- c(
  "age", 
  "sexe", 
  "scoreKnaus", 
  "glasgowIoa", 
  "crpAdmission", 
  "creatinineAdmission", 
  "leucocytesTotadmission", 
  "dfg", 
  "obesiteSurpoids", 
  "presentationCliniqueUrgences"
)

df_ebal <- df_ipw %>%
  select(A, Y, all_of(vars_ebal))

w.ebal <- weightit(A ~ ., data = df_ebal, method = "ebal", estimand = "ATE")

df_ebal$poids_ebal <- get.w(w.ebal)

bal_ebal <- bal.tab(w.ebal, un = TRUE, m.threshold = 0.1)
print(bal_ebal)

#Modèle de régression pondéré (EB)
modele_ebal <- glm(Y ~ A, data = df_ebal, weights = poids_ebal, family = binomial)

effet_ebal <- tidy(modele_ebal, exponentiate = TRUE, conf.int = TRUE)

effet_ebal_A <- effet_ebal %>%
  filter(term == "A") %>%
  mutate(
    Variable = "A (antibio)",
    OR = round(estimate, 2),
    `IC 95%` = paste0("[", round(conf.low, 2), " ; ", round(conf.high, 2), "]"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  select(Variable, OR, `IC 95%`, `p-value`)

print(effet_ebal_A)


###############################################################
##############################################################
##GRAPHIQUE:

#1. Tableau comparatif pondéré vs non pondéré
tab_compare <- bal.tab(w.out, un = TRUE, disp.v.ratio = TRUE)

# Convertir en data frame
balance_df <- as.data.frame(tab_compare$Balance)
View(balance_df)  


# Extraire les coefficients exponentiés et leurs IC
effet <- tidy(modele_iptw, exponentiate = TRUE, conf.int = TRUE)

# Filtrer uniquement la variable de traitement "A"
effet_A <- effet %>% filter(term == "A")

# Créer le forest plot
ggplot(effet_A, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#0072B2") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    x = "Odds Ratio (IC 95%)",
    y = "",
    title = "Effet marginal du traitement (IPTW)",
    subtitle = paste0("OR = ", round(effet_A$estimate, 3), 
                      " [", round(effet_A$conf.low, 3), "; ", round(effet_A$conf.high, 3), "]")
  ) +
  theme_minimal(base_size = 14)


#######################################################################
############TABLEAU COMPARATIF
# Analyse brute
effet_brut_A <- tidy(modele_brut, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "A") %>%
  mutate(
    Méthode = "Analyse brute",
    OR = round(estimate, 2),
    `IC 95%` = paste0("[", round(conf.low, 2), " ; ", round(conf.high, 2), "]"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  mutate(`p-value` = as.character(`p-value`)) %>%  # uniformiser le type
  select(Méthode, OR, `IC 95%`, `p-value`)

# IPTW
effet_iptw_A <- tidy(modele_iptw, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "A") %>%
  mutate(
    Méthode = "IPTW",
    OR = round(estimate, 2),
    `IC 95%` = paste0("[", round(conf.low, 2), " ; ", round(conf.high, 2), "]"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  mutate(`p-value` = as.character(`p-value`)) %>%
  select(Méthode, OR, `IC 95%`, `p-value`)

# Entropy Balancing
effet_ebal_A <- tidy(modele_ebal, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "A") %>%
  mutate(
    Méthode = "Entropy Balancing",
    OR = round(estimate, 2),
    `IC 95%` = paste0("[", round(conf.low, 2), " ; ", round(conf.high, 2), "]"),
    `p-value` = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))
  ) %>%
  mutate(`p-value` = as.character(`p-value`)) %>%
  select(Méthode, OR, `IC 95%`, `p-value`)

# Fusionner les tableaux
tableau_comparatif <- bind_rows(effet_brut_A, effet_iptw_A, effet_ebal_A)

# Affichage
print(tableau_comparatif)
View(tableau_comparatif)


# Recréer les données pour le graphique avec les bornes IC numériques
effets_plot <- bind_rows(
  tidy(modele_brut, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Méthode = "Analyse brute"),
  tidy(modele_iptw, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Méthode = "IPTW"),
  tidy(modele_ebal, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "A") %>% mutate(Méthode = "Entropy Balancing")
)

ggplot(effets_plot, aes(x = estimate, y = Méthode)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "#0072B2") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "Effet du traitement selon la méthode",
    x = "Odds Ratio (IC 95%)",
    y = "Méthode"
  ) +
  theme_minimal(base_size = 14)

