# Statistical Analysis of Simpson's Paradox and Log-Linear Models in R

This repository contains the implementation and comprehensive statistical analysis of a classic problem in **Categorical Data Analysis** using **R**. The primary objective of this project is to investigate **Simpson's Paradox** in simulated clinical data (comparing the efficacy of two cholesterol-lowering drugs) using **marginal and conditional contingency tables**, **Odds Ratios (OR)**, and **Log-Linear Models**.

## 📋 Table of Contents

1.  [Problem Overview & Data Structure](#1-problem-overview--data-structure)
2.  [Executed R Code](#2-executed-r-code)
3.  [Output Analysis & Answers to Five Key Questions](#3-output-analysis--answers-to-five-key-questions)
    -   [Question 1: Marginal Table — Which Drug Is Better Overall?](#question-1-marginal-table--which-drug-is-better-overall)
    -   [Question 2: Conditional Tables — Stratified by Cholesterol Level](#question-2-conditional-tables--stratified-by-cholesterol-level)
    -   [Question 3: Did Simpson's Paradox Occur?](#question-3-did-simpsons-paradox-occur)
    -   [Question 4: Log-Linear Model Fitting (`loglm`) & Odds Ratio Interpretation](#question-4-log-linear-model-fitting-loglm--odds-ratio-interpretation)
    -   [Question 5: ANOVA & Evaluation of the Three-Way Interaction](#question-5-anova--evaluation-of-the-three-way-interaction)
4.  [Final Conclusions](#4-final-conclusions)

------------------------------------------------------------------------

## 📊 1. Problem Overview & Data Structure

In this study, $N = 500$ patients were treated with either **Drug A** or **Drug B**. The three categorical variables are:

-   **Treatment (**$T$): Two levels — `Drug A` ($n = 250$) and `Drug B` ($n = 250$).
-   **Baseline Cholesterol Level (**$C$ - `Chol_Level`): Confounding variable with two levels — `High` and `Normal`.
-   **Treatment Outcome (**$O$ - `Outcome`): Two levels — `Reduced` (success/cholesterol reduction) and `Not Reduced` (failure).

### Summary of the Three-Way Contingency Table (`ftable`)

| Treatment | Cholesterol Level (`Chol_Level`) | `Reduced` | `Not Reduced` | Subgroup Size ($n$) | Success Rate |
|:-----------|:-----------|:----------:|:----------:|:----------:|:----------:|
| **Drug A** | High | 160 | 40 | 200 | **80.0%** |
| **Drug A** | Normal | 30 | 20 | 50 | **60.0%** |
| **Drug B** | High | 45 | 5 | 50 | **90.0%** |
| **Drug B** | Normal | 130 | 70 | 200 | **65.0%** |

------------------------------------------------------------------------

## 💻 2. Executed R Code

The `MASS` package is required to fit the hierarchical log-linear models via `loglm()`:

``` r
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
                 rep("Reduced", 130), rep("Not Reduced", 70))   # B | Normal : 65%
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

# Mutual independence: [T][C][O]
l1 <- loglm(~ Treatment + Chol_Level + Outcome, data = med_table)

# Conditional independence: [TC][CO]
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
```

------------------------------------------------------------------------

## 🔍 3. Output Analysis & Answers to Five Key Questions

### Question 1: Marginal Table — Which Drug Is Better Overall?

Ignoring baseline cholesterol levels and collapsing the data into a $2 \times 2$ marginal table ($T \times O$):

| Treatment  | `Reduced` | `Not Reduced` | Total | Success Rate |
|:-----------|:---------:|:-------------:|:-----:|:------------:|
| **Drug A** |    190    |      60       |  250  |  **76.0%**   |
| **Drug B** |    175    |      75       |  250  |  **70.0%**   |

-   **Marginal Odds Ratio (**$\text{OR}_{\text{marginal}}$):

$$\text{OR} = \frac{190 \times 75}{60 \times 175} = \mathbf{1.357}$$

*(The apparent odds of cholesterol reduction with Drug A are roughly* $1.36$ times those of Drug B).

-   **Significance Testing:** Both the Pearson Chi-Square test ($\chi^2 \approx 2.28, \; df = 1, \; p \approx 0.13$) and the Likelihood Ratio test ($\Delta G^2 = 2.40, \; df = 1, \; p = 0.121$) indicate that this difference is **not statistically significant**.

> ✅ **Answer:** Overall (in the marginal table), **Drug A** appears to perform better ($76\%$ vs. $70\%$), though this apparent advantage is not statistically significant.

------------------------------------------------------------------------

### Question 2: Conditional Tables — Stratified by Cholesterol Level

Controlling for baseline cholesterol (`Chol_Level`) and examining each stratum separately:

#### A) High Cholesterol Stratum (`High`):

| Treatment  | `Reduced` | `Not Reduced` | Success Rate |
|:-----------|:---------:|:-------------:|:------------:|
| **Drug A** |    160    |      40       |    80.0%     |
| **Drug B** |    45     |       5       |  **90.0%**   |

-   $\text{OR}_{A \text{ vs } B \mid \text{High}} = \frac{160 \times 5}{40 \times 45} = \mathbf{0.444}$ (equivalently, $\text{OR}_{B \text{ vs } A} = \mathbf{2.25}$). **Drug B** outperforms Drug A in this subgroup ($p \approx 0.11$).

#### B) Normal Cholesterol Stratum (`Normal`):

| Treatment  | `Reduced` | `Not Reduced` | Success Rate |
|:-----------|:---------:|:-------------:|:------------:|
| **Drug A** |    30     |      20       |    60.0%     |
| **Drug B** |    130    |      70       |  **65.0%**   |

-   $\text{OR}_{A \text{ vs } B \mid \text{Normal}} = \frac{30 \times 70}{20 \times 130} = \mathbf{0.808}$ (equivalently, $\text{OR}_{B \text{ vs } A} = \mathbf{1.24}$). **Drug B** also outperforms Drug A in this subgroup ($p \approx 0.51$).

-   **Common Mantel–Haenszel Odds Ratio (**$\text{OR}_{\text{MH}}$): Equal to $\mathbf{0.66}$ for A vs. B (or $\mathbf{1.51}$ favoring Drug B over Drug A across both strata).

> ✅ **Answer:** After controlling for baseline cholesterol level, **Drug B outperforms Drug A in both subgroups.**

------------------------------------------------------------------------

### Question 3: Did Simpson's Paradox Occur? {#question-3-did-simpsons-paradox-occur}

> ✅ **Yes — Simpson's Paradox occurred descriptively.**

-   **Marginally:** $\text{Drug A } (76\%) > \text{Drug B } (70\%)$
-   **Conditionally:** Within both cholesterol strata, $\text{Drug B} > \text{Drug A}$ ($90\% > 80\%$ in `High`, and $65\% > 60\%$ in `Normal`)

#### Mechanism & Cause (Confounding):

1.  **Severe Imbalance in Treatment Allocation (**$\text{OR}_{T \times C} = 16$): $80\%$ of patients receiving Drug A had high cholesterol ($200/250$), whereas only $20\%$ of patients receiving Drug B had high cholesterol ($50/250$).
2.  **Strong Prognostic Effect of Cholesterol Level (**$\text{OR}_{C \times O} = 2.56, \; p < 0.001$): Patients with high baseline cholesterol experienced a much higher overall reduction rate ($82\%$) compared to patients with normal cholesterol ($64\%$).
3.  **Weighted Average Effect:**

$$\text{Success}_A = (0.80 \times 80\%) + (0.20 \times 60\%) = \mathbf{76.0\%}$$

$$\text{Success}_B = (0.20 \times 90\%) + (0.80 \times 65\%) = \mathbf{70.0\%}$$

Thus, the marginal superiority of Drug A is entirely a **confounding artifact** rather than a true pharmacological advantage.

------------------------------------------------------------------------

### Question 4: Log-Linear Model Fitting (`loglm`) & Odds Ratio Interpretation

Four hierarchical log-linear models were fitted to the 3-way contingency table:

| Model | Structure | Formula | Deviance ($G^2$) | $df$ | $P$-value | Goodness of Fit |
|:---------:|:----------|:---------:|:---------:|:---------:|:---------:|:----------|
| **l1** | Mutual Independence | `[T][C][O]` | $217.055$ | $4$ | $< 0.001$ | ❌ Rejected |
| **l2** | Conditional Independence | `[TC][CO]` | $3.457$ | $2$ | $0.178$ | ✅ **Excellent fit (Optimal)** |
| **l3** | Homogeneous Association | `[TC][CO][TO]` | $1.054$ | $1$ | $0.305$ | ✅ Good fit |
| **l4** | Saturated | `[TCO]` | $0.000$ | $0$ | $1.000$ | Exact fit |

-   **Interpretation of Selected Model (`l2`):** Model `l2` represents **conditional independence between Treatment and Outcome given Cholesterol Level** ($T \perp O \mid C$). Because `l2` fits the data well ($p = 0.178$), it confirms that once cholesterol level is accounted for, no statistically significant association remains between drug choice and outcome. The combination of *marginal association* and *conditional independence* is the exact log-linear formulation of Simpson's Paradox.

#### Computing True Odds Ratios from `exp(coef(l3))`:

The `loglm()` function in R uses **sum-to-zero constraints** ($\lambda_{11} = -\lambda_{12} = -\lambda_{21} = \lambda_{22}$). Consequently, the values returned by `exp(coef(l3))` represent the **fourth root of the Odds Ratio** and must be raised to the **4th power** to obtain the model-based ORs:

1.  **Treatment** $\times$ Cholesterol (`Treatment:Chol_Level`):

$$\text{OR}_{TC} = (2.046)^4 \approx \mathbf{17.52}$$

*(Confirms the heavy allocation of high-cholesterol patients to Drug A).*

2.  **Cholesterol** $\times$ Outcome (`Chol_Level:Outcome`):

$$\text{OR}_{CO} = (1.345)^4 \approx \mathbf{3.27}$$

*(The odds of cholesterol reduction are roughly* $3.27$ times higher in the High cholesterol group than in the Normal group).

3.  **Conditional Treatment** $\times$ Outcome (`Treatment:Outcome`):

$$\text{OR}_{B \text{ vs } A \mid C} = (1.106)^4 \approx \mathbf{1.50}$$

$$\text{OR}_{A \text{ vs } B \mid C} = (0.904)^4 \approx \mathbf{0.67}$$

*(Reflects the true direction of the effect: controlling for cholesterol, the odds of reduction with Drug B are roughly* $1.50$ times those of Drug A).

------------------------------------------------------------------------

### Question 5: ANOVA & Evaluation of the Three-Way Interaction

Comparison of nested models using Likelihood Ratio Tests (`anova(l1, l2, l3, l4)`):

| Model Comparison | $\Delta G^2$ (`Delta Dev`) | $\Delta df$ | $P$-value | Hypothesis Tested |
|:-------------:|:-------------:|:-------------:|:-------------:|:--------------|
| $l_1 \rightarrow l_2$ | $213.598$ | $2$ | $< 0.001$ | Necessity of confounding associations `[TC]` and `[CO]` |
| $l_2 \rightarrow l_3$ | $2.402$ | $1$ | $0.121$ | Significance of direct `[TO]` association controlling for cholesterol |
| $l_3 \rightarrow l_4$ | $1.054$ | $1$ | $0.305$ | Significance of the **three-way interaction `[TCO]`** |

#### ⚠️ Clarifying a Common Misconception About the Three-Way Interaction (`TCO`):

Contrary to the assumption that proving Simpson's Paradox requires a statistically significant three-way interaction (`TCO`), the ANOVA output demonstrates:

> **The three-way interaction (`TCO`) is NOT statistically significant (**$\Delta G^2 = 1.054, \; df = 1, \; p = 0.305$).

-   **Statistical Meaning:** The **homogeneous association** assumption holds—meaning the conditional odds ratios comparing Drug B to Drug A in the High cholesterol group ($\text{OR} = 2.25$) and Normal cholesterol group ($\text{OR} = 1.24$) do not significantly differ from one another (particularly given the small sample sizes of $n = 50$ in the `A-Normal` and `B-High` cells, which limit statistical power).

-   **How Simpson's Paradox Is Demonstrated Instead:** Simpson's Paradox does not require a three-way interaction; it is established through **three complementary pieces of evidence**:

    1.  **Descriptive Direction Reversal:** Marginal superiority of Drug A ($76\% > 70\%$) vs. conditional superiority of Drug B in every stratum ($90\% > 80\%$ and $65\% > 60\%$).
    2.  **Adequacy of the Conditional Independence Model (`l2`):** Model `[TC][CO]` fits well ($p = 0.178$), showing that marginal association disappears once the confounder is conditioned on.
    3.  **Highly Significant Confounding Pathways:** Adding `[TC]` and `[CO]` reduces deviance by $\Delta G^2 = 213.6$ ($p < 0.001$).

------------------------------------------------------------------------

## 🎯 4. Final Conclusions

1.  **Most Parsimonious Model:** Model **`l2` (`~ Treatment * Chol_Level + Chol_Level * Outcome`)** is the simplest and best-fitting model for these data.
2.  **Nature of the Paradox:** The dataset exhibits a **descriptive Simpson's Paradox driven by strong confounding**. Drug B is superior within both patient strata, but unbalanced treatment assignment of high-cholesterol patients to Drug A reverses the marginal comparison.
3.  **Statistical Recommendation:** To establish statistical significance for the conditional superiority of Drug B (`l3`) or stratum-specific effect modification (`l4`), larger sample sizes are needed in the underrepresented cells (`Drug A | Normal` and `Drug B | High`).
