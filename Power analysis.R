## Power analysis

# Load necessary libraries
#install.packages("simstudy")

# Load necessary libraries
library(simstudy)
library(broom)

# Set up parameters
n_simulations <- 1000   # Number of simulations
n_sample <- 1380         # Sample size per simulation, adjust as needed
p_outcome <- 0.5       # Overall probability of the outcome P(Y=1) (dementia 0.15, cognitive impairment 0.5)
p_exposure <- 0.35      # Probability of the exposure P(X=1)
odds_ratio <- 1.4         # Assumed odds ratio for the exposure effect
alpha <- 0.05           # Significance level

# Calculate log-odds ratio
log_odds_ratio <- log(odds_ratio)

# Calculate baseline log-odds for y when x=0
p_y_given_x0 <- p_outcome / (1 - p_exposure * (odds_ratio - 1))  # Adjusted for overall probability
log_odds_y_given_x0 <- log(p_y_given_x0 / (1 - p_y_given_x0))

# Define the simulation function incorporating outcome probability
simulate_data <- function(n_sample, p_exposure, log_odds_y_given_x0, log_odds_ratio, alpha) {
  
  # Define the data for exposure
  def <- defData(varname = "x", dist = "binary", formula = p_exposure)
  
  # Expand data based on defined exposure probabilities
  data <- genData(n_sample, def)
  
  # Calculate the logit for the outcome with exposure effect
  data <- data[, logit_p := log_odds_y_given_x0 + log_odds_ratio * x]
  
  # Convert logit to probability
  data <- data[, p := exp(logit_p) / (1 + exp(logit_p))]
  
  # Simulate the outcome y based on calculated probabilities
  data <- data[, y := rbinom(n_sample, 1, p)]
  
  # Fit logistic regression model
  model <- glm(y ~ x, data = data, family = binomial())
  
  # Extract p-value for the exposure effect
  p_value <- summary(model)$coefficients["x", "Pr(>|z|)"]
  
  return(p_value)
}

# Perform the simulations
set.seed(123)  # For reproducibility
p_values <- replicate(n_simulations, simulate_data(n_sample, p_exposure, log_odds_y_given_x0, log_odds_ratio, alpha))

# Calculate observed power
observed_power <- mean(p_values < alpha)
print(paste("Observed Power:", observed_power))




