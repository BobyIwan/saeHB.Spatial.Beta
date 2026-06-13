## code to prepare `DATASET` dataset goes here

set.seed(20260601L)

# Setup dimensions and main parameters

n_row <- 6L
n_col <- 6L
m     <- n_row * n_col

beta_true <- c(`(Intercept)` = 1, x1 = 1, x2 = 1)
rho_true  <- 0.70
sigma_u2  <- 1.0

domain_names <- sprintf("area_%03d", seq_len(m))
eps <- 1e-6

# Generate spatial matrices
nb_queen <- spdep::cell2nb(
  n_row,
  n_col,
  type = "queen"
)

# Row-standardized spatial weight matrix (SAR and Moran's I)
weight_mat <- spdep::nb2mat(
  nb_queen,
  style = "W",
  zero.policy = TRUE
)

# Binary adjacency matrix (CAR Leroux)
adjacency_mat <- spdep::nb2mat(
  nb_queen,
  style = "B",
  zero.policy = TRUE
)

rownames(weight_mat) <- colnames(weight_mat) <- domain_names
rownames(adjacency_mat) <- colnames(adjacency_mat) <- domain_names

# Generate auxiliary variables
x1 <- stats::rnorm(m)
x2 <- stats::rnorm(m)

X <- cbind(1, x1, x2)

# Generate design-effect information
n_i <- round(stats::runif(m, 10, 50))
deff <- round(stats::runif(m, 1, 2.5), 2)

phi_deff <- (n_i / deff) - 1

# Generate SAR spatial random effects
u <- stats::rnorm(
  m,
  mean = 0,
  sd = sqrt(sigma_u2)
)

v <- as.numeric(
  solve(diag(m) - rho_true * weight_mat) %*% u
)

# Generate mean parameter
mu <- stats::plogis(
  as.numeric(X %*% beta_true + v)
)

# Generate beta response
y <- pmin(
  pmax(
    stats::rbeta(
      m,
      mu * phi_deff,
      (1 - mu) * phi_deff
    ),
    eps
  ),
  1 - eps
)

# Complete dataset
databeta <- data.frame(
  domain = domain_names,
  y = y,
  x1 = x1,
  x2 = x2,
  n_i = n_i,
  deff = deff,
  stringsAsFactors = FALSE
)

# Dataset with missing responses
databeta_na <- databeta

databeta_na$y[c(7, 13, 18, 21, 32)] <- NA

# Export datasets
usethis::use_data(
  databeta,
  databeta_na,
  weight_mat,
  adjacency_mat,
  overwrite = TRUE,
  compress = "xz"
)
