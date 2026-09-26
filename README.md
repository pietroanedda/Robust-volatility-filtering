# NASDAQ-100 Volatility Analysis (2020-2024)

This repository contains the code and the final report for an advanced time series analysis project. The main goal is to study the variability and the volatility dynamics of the NASDAQ-100 index during the 2020-2024 period.

## Project Description
The analysis focuses on the daily log-returns of the NASDAQ-100. Empirical data highlights the "volatility clustering" phenomenon and heavy leptokurtosis (heavy tails), which invalidate the classic assumption of normally distributed returns. 

To model these dynamics, three different models were implemented, estimated, and compared:
* **Gaussian GARCH**: The traditional model used as a baseline for comparison.
* **Beta-t-GARCH**: A model belonging to the newer class of score-driven models, based on a conditional Student's t-distribution.
* **Beta-t-EGARCH**: A second score-driven model that models the logarithm of the variance, ensuring its positivity without strict constraints on the static parameters.

## Repository Structure
The repository includes the R scripts used to conduct the estimations and diagnostic tests:

* `GARCH.R`: Contains functions to compute the log-likelihood, filter the dynamic scale parameter, and estimate the Gaussian GARCH model. It implements a penalized likelihood approach to enforce stationarity constraints using the `nlminb` optimizer.
* `beta_t_GARCH.R`: Script containing functions to implement the Beta-t-GARCH model. It includes the definition of the conditional score, the likelihood function, and the estimation of optimal parameters (including the degrees of freedom for the Student's t-distribution).
* `beta-t-EGARCH.R`: Contains functions to estimate, filter, and calculate the standardized residuals for the Beta-t-EGARCH model.
* `stat_functions.R`: Script containing auxiliary statistical functions. Specifically, it includes `acf_test` for Ljung-Box/Box-Pierce tests, `qqplot_t` for plotting theoretical Q-Q plots for the Student's t-distribution, and `IC` for computing the AIC and BIC.

## Data Source
* Data was collected using the `quantmod` R package, sourcing directly from Yahoo Finance.
* The analysis was conducted on daily adjusted close prices to account for dividend payouts and stock splits.

## Main Results
* The Gaussian GARCH model proved inadequate in capturing extreme market drawdowns (e.g., the COVID-19 pandemic and the outbreak of the Russo-Ukrainian war) due to its flawed distributional assumptions.
* The score-driven models proved to be more robust filters against extreme values, as their innovations (based on conditional scores) are uniformly bounded.
* Model selection (via AIC and BIC criteria) alongside diagnostic checks identified the **Beta-t-EGARCH** as the best model among the three. It successfully ensures that both standardized residuals and score innovations behave as martingale difference sequences, validating its theoretical optimality under heavy-tailed data.