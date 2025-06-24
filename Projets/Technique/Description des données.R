install.packages("medicaldata")
library(medicaldata)
install.packages("psych")
library(psych)
install.packages("naniar")
library(naniar)
BDD <- cytomegalovirus
str(BDD)
BDD$sex <- factor(BDD$sex, labels = c ("Female" , "Male"))
BDD$sex
str(BDD)
BDD$race <- factor(BDD$race, labels = c("African American", "White"))
BDD$prior.radiation <- as.logical(BDD$prior.radiation, labels= c ("No", "Yes"))
BDD$prior.transplant <- as.logical(BDD$prior.transplant, labels= c ("No", "Yes"))
BDD$donor.sex <- factor(BDD$donor.sex, labels = c ("Female" , "Male"))
BDD$cmv <- as.logical(BDD$cmv, labels= c ("No", "Yes"))
BDD$agvhd <- as.logical(BDD$agvhd, labels= c ("No", "Yes"))
BDD$cgvhd <- as.logical(BDD$cgvhd, labels= c ("No", "Yes"))
BDD$diagnosis.type <- factor(BDD$diagnosis.type, labels = c("Lymphoid", "Myeloid"))
BDD$`C1/C2`<- factor(BDD$`C1/C2`, labels = c("Heterozygous", "Homozygous"))
BDD$recipient.cmv <- as.logical(BDD$recipient.cmv, labels= c ("Negative", "Positive"))
BDD$donor.cmv <- as.logical(BDD$donor.cmv, labels= c ("Negative", "Positive"))
str(BDD)
describe(BDD)

describe_BDD <- function(BDD){
  for (col_name in names(BDD)) {
    cat("Variable:", col_name, "\n")
    column <- BDD[[col_name]]
    
    if (is.numeric(column)) {
      cat("Type: Quantitative\n")
      cat(" Moyenne:", mean(column, na.rm = TRUE), "\n")
      cat(" Médiane:", median(column, na.rm = TRUE), "\n")
      cat(" Écart-type:", sd(column, na.rm = TRUE), "\n")
      cat(" Étendue:", range(column, na.rm = TRUE), "\n") 
      cat(" IQR:", IQR(column, na.rm = TRUE), "\n\n")
      
    } 
    
    else if (is.factor(column) || is.character(column) || is.logical(column)) {
      cat("Type: Qualitative\n")
      cat("  Fréquences:\n")
      print(table(column, useNA = "ifany"))
      cat("\n")
      
    } 
    else {
      cat("Type: Non pris en charge\n\n")
    }
  }
}
      
describe_BDD(BDD)
      