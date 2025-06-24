library(ggplot2)
library(readxl)
library(medicaldata)
library(dplyr)
library(tidyr)
data <- cytomegalovirus

#Modification des variables
data$sex <- factor(data$sex, labels = c ("Female" , "Male"))
data$diagnosis.type <- factor(data$diagnosis.type, labels = c("Lymphoid", "Myeloid"))
data$race <- factor(data$race, labels = c("African American", "White"))
data_filtered <- data %>% filter(!is.na(diagnosis.type))

#Graphe 1: Boxplot: distribution de l'age selon type de diagnostic et sexe (comparaison)
ggplot(data_filtered, aes(x = diagnosis.type, y = age, fill = sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("pink", "skyblue")) + # Palette de couleurs
  labs(
    title = "Distribution de l'âge des patients selon le type de diagnostic et le sexe",
    x = "Type de Diagnostic",
    y = "Âge des patients",
    fill = "Sex"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Graphe 2: nuages de points
ggplot(data, aes(x = age, y = time.to.transplant, color = sex)) +
  geom_point() +
  facet_wrap(~ race) +
  scale_color_manual(values = c("Female" = "deeppink", 
                                "Male" = "blue")) + 
  labs(title = "Temps de transplantation en fonction de l'age, race et sexe",
       x = "Age", y = "Time to Transplant") +
  theme_minimal()
