set.seed(123)

# Simpson's Paradox demonstration data
# Drug B is better within each cholesterol stratum, but appears worse marginally
med_data <- data.frame(
  Treatment  = c(rep("Drug A", 200), rep("Drug A", 50),
                 rep("Drug B", 50),  rep("Drug B", 200)),
  Chol_Level = c(rep("High", 200),   rep("Normal", 50),
                 rep("High", 50),    rep("Normal", 200)),
  Outcome    = c(rep("Reduced", 160), rep("Not Reduced", 40),   # A | High   : 80%
                 rep("Reduced", 30),  rep("Not Reduced", 20),   # A | Normal : 60%
                 rep("Reduced", 45),  rep("Not Reduced", 5),    # B | High   : 90%
                 rep("Reduced", 130), rep("Not Reduced", 70))    # B | Normal : 65%
)

# 3-way contingency table
med_table <- table(Treatment  = med_data$Treatment,
                   Chol_Level = med_data$Chol_Level,
                   Outcome    = med_data$Outcome)

ftable(med_table)


# ---------- 1) Marginal Analysis ----------
marg_table <- margin.table(med_table, margin = c(1, 3))

cat("\n=== Marginal Table (Treatment x Outcome) ===\n")
print(marg_table)

cat("\n=== Row Proportions (by Treatment) ===\n")
print(prop.table(marg_table, margin = 1))

cat("\n=== Marginal Success Rate (Reduced) ===\n")
print(round(prop.table(marg_table, margin = 1)[, "Reduced"], 3))


# ---------- 2) Conditional Analysis (given Chol_Level) ----------
cat("\n=== Conditional Tables (Treatment x Outcome | Chol_Level) ===\n")

for (lev in dimnames(med_table)$Chol_Level) {
  cat("\n--- Chol_Level =", lev, "---\n")
  cond_tab <- med_table[, lev, ]
  print(cond_tab)
  
  cat("Conditional Success Rate (Reduced):\n")
  print(round(prop.table(cond_tab, margin = 1)[, "Reduced"], 3))
}


# ---------- 3) Log-Linear Models ----------
library(MASS)

# Mutual independence
l1 <- loglm(~ Treatment + Chol_Level + Outcome, data = med_table)

# Conditional independence
l2 <- loglm(~ Treatment * Chol_Level + Chol_Level * Outcome, data = med_table)

# Homogeneous association: [TC][CO][TO]
l3 <- loglm(~ (Treatment + Chol_Level + Outcome)^2, data = med_table)

# Saturated: [TCO]
l4 <- loglm(~ Treatment * Chol_Level * Outcome, data = med_table)


# ---------- 4) Model Comparison ----------
cat("\n=== Model Comparison (ANOVA) ===\n")
print(anova(l1, l2, l3, l4))


# ---------- 5) Odds Ratios from Homogeneous Association Model ----------
cat("\n=== Odds Ratios (Model l3) ===\n")

cat("\nTreatment x Chol_Level:\n")
print(round(exp(coef(l3)$Treatment.Chol_Level), 3))

cat("\nChol_Level x Outcome:\n")
print(round(exp(coef(l3)$Chol_Level.Outcome), 3))

cat("\nTreatment x Outcome:\n")
print(round(exp(coef(l3)$Treatment.Outcome), 3))

