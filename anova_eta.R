anova_eta <- function(model, round_by) {
  library(effectsize)
  
  a <- anova(model)
  b <- eta_squared(model)
  
  sign <- c()
  for (i in 1:nrow(a)) {
  if (a$`Pr(>F)`[i] < 0.001) {
    sign <- append(sign, "***")}
    else if (a$`Pr(>F)`[i] < 0.01) {
      sign <- append(sign, "**") }
    else if (a$`Pr(>F)`[i] < 0.05) {
      sign <- append(sign, "*") }
      else {sign <- append(sign, "")}}
  
  c <- data.frame(round(a,round_by), sign, round(b$Eta2_partial, round_by), 
                  paste("[", round(b$CI_low, round_by), ", ", round(b$CI_high, round_by), "]"))
  colnames(c) <- c("SS", "MS", "N.DF", "D.DF", "F", "p", "", "Eta2 par.", "95%CI")
  
  return(c)
}