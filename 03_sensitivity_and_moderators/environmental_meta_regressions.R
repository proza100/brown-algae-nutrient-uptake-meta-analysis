options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Usage: Rscript environmental_meta_regressions.R effects.csv temperature_context.csv irradiance_context.csv output_directory")
input_effects <- args[1]
temp_context <- args[2]
irr_context <- args[3]
out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

effects <- read.csv(input_effects, check.names = FALSE)
temp <- read.csv(temp_context, check.names = FALSE)
irr <- read.csv(irr_context, check.names = FALSE)
names(temp) <- sub("^\\ufeff", "", names(temp))
names(irr) <- sub("^\\ufeff", "", names(irr))
if (grepl("Study_ID$", names(temp)[1])) names(temp)[1] <- "Study_ID"
if (grepl("Study_ID$", names(irr)[1])) names(irr)[1] <- "Study_ID"

needed <- c("Effect_ID", "Study_ID", "Experiment_ID", "Nutrient_form", "y", "v")
stopifnot(all(needed %in% names(effects)))

prepare_context <- function(context, driver) {
  keep <- effects[effects$Driver_harmonised == driver, needed]
  merged <- merge(context, keep, by = c("Effect_ID", "Study_ID", "Experiment_ID"), all.x = TRUE)
  if (any(is.na(merged$y)) || any(is.na(merged$v))) {
    stop(paste("Unmatched analytical effects for", driver))
  }
  merged
}

temp <- prepare_context(temp, "Temperature")
irr <- prepare_context(irr, "Irradiance")

covariance_block <- function(group, rho, tau_study, tau_experiment) {
  v <- as.numeric(group$v)
  exp_id <- as.character(group$Experiment_ID)
  same_exp <- outer(exp_id, exp_id, "==") * 1
  sampling <- diag(v, nrow = length(v), ncol = length(v)) +
    rho * sqrt(outer(v, v)) * (same_exp - diag(nrow(group)))
  sampling + tau_study * matrix(1, nrow(group), nrow(group)) + tau_experiment * same_exp
}

gls_components <- function(data, rho, tau_study, tau_experiment) {
  p <- 2L
  xtwx <- matrix(0, p, p)
  xtwy <- numeric(p)
  logdet <- 0
  quad <- 0
  pieces <- list()
  split_data <- split(data, data$Study_ID)
  for (study_id in names(split_data)) {
    g <- split_data[[study_id]]
    y <- as.numeric(g$y)
    X <- cbind(1, as.numeric(g$moderator_scaled_centered))
    V <- covariance_block(g, rho, tau_study, tau_experiment)
    R <- try(chol(V), silent = TRUE)
    if (inherits(R, "try-error")) return(NULL)
    W <- chol2inv(R)
    xtwx <- xtwx + crossprod(X, W %*% X)
    xtwy <- xtwy + crossprod(X, W %*% y)[, 1]
    logdet <- logdet + 2 * sum(log(diag(R)))
    quad <- quad + as.numeric(crossprod(y, W %*% y))
    pieces[[length(pieces) + 1L]] <- list(group = g, X = X, W = W)
  }
  beta <- try(solve(xtwx, xtwy), silent = TRUE)
  if (inherits(beta, "try-error")) return(NULL)
  det_x <- determinant(xtwx, logarithm = TRUE)
  if (det_x$sign <= 0) return(NULL)
  reml <- logdet + as.numeric(det_x$modulus) + quad - as.numeric(crossprod(beta, xtwy))
  list(beta = beta, xtwx = xtwx, reml = reml, pieces = pieces)
}

fit_model <- function(data, rho) {
  k <- length(unique(data$Study_ID))
  p <- 2L
  if (k <= p) stop("Too few studies")
  scale0 <- max(var(data$y), median(data$v), 0.01, na.rm = TRUE)
  objective <- function(par) {
    ans <- gls_components(data, rho, par[1], par[2])
    if (is.null(ans) || !is.finite(ans$reml)) return(1e100)
    0.5 * ans$reml
  }
  starts <- rbind(c(0, 0), c(scale0 / 4, scale0 / 4), c(scale0, 0), c(0, scale0), c(scale0, scale0))
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    optim(starts[i, ], objective, method = "L-BFGS-B", lower = c(0, 0),
          control = list(maxit = 1000, factr = 1e7))
  })
  best <- fits[[which.min(vapply(fits, function(x) x$value, numeric(1)))]]
  tau <- pmax(best$par, 0)
  ans <- gls_components(data, rho, tau[1], tau[2])
  beta <- ans$beta
  bread <- solve(ans$xtwx)
  meat <- matrix(0, p, p)
  for (piece in ans$pieces) {
    residual <- as.numeric(piece$group$y) - as.numeric(piece$X %*% beta)
    score <- crossprod(piece$X, piece$W %*% residual)[, 1]
    meat <- meat + tcrossprod(score)
  }
  robust <- (k / (k - p)) * bread %*% meat %*% bread
  se_model <- sqrt(pmax(diag(bread), 0))
  se_cluster <- sqrt(pmax(diag(robust), 0))
  se_used <- pmax(se_model, se_cluster)
  df <- k - p
  crit <- qt(0.975, df)
  lo <- beta - crit * se_used
  hi <- beta + crit * se_used
  p_value <- 2 * pt(-abs(beta / se_used), df)
  list(beta = beta, se_model = se_model, se_cluster = se_cluster, se_used = se_used,
       lo = lo, hi = hi, p_value = p_value, tau = tau, df = df,
       convergence = best$convergence == 0)
}

specifications <- data.frame(
  Driver = c(rep("Temperature", 3), rep("Irradiance", 2)),
  Moderator = c("Higher_temperature_C", "Temperature_contrast_C", "Higher_minus_acclimation_C",
                "Higher_irradiance_PPFD", "Irradiance_contrast_PPFD"),
  Scale = c(1, 1, 1, 100, 100),
  Unit = c("per 1 degree C", "per 1 degree C", "per 1 degree C",
           "per 100 umol photons m-2 s-1", "per 100 umol photons m-2 s-1")
)

result_rows <- list()
feasibility_rows <- list()
loo_rows <- list()

for (s in seq_len(nrow(specifications))) {
  spec <- specifications[s, ]
  source <- if (spec$Driver == "Temperature") temp else irr
  for (nutrient in intersect(c("Nitrate", "Ammonium", "Phosphate"), unique(source$Nutrient))) {
    g <- source[source$Nutrient == nutrient, ]
    complete <- g[is.finite(g[[spec$Moderator]]), ]
    studies <- length(unique(complete$Study_ID))
    experiments <- length(unique(complete$Experiment_ID))
    feasibility_rows[[length(feasibility_rows) + 1L]] <- data.frame(
      Driver = spec$Driver, Nutrient = nutrient, Moderator = spec$Moderator,
      effects = nrow(complete), studies = studies, experiments = experiments,
      min = if (nrow(complete)) min(complete[[spec$Moderator]]) else NA,
      median = if (nrow(complete)) median(complete[[spec$Moderator]]) else NA,
      max = if (nrow(complete)) max(complete[[spec$Moderator]]) else NA,
      feasible = studies >= 3
    )
    if (studies < 3) next
    scaled <- complete[[spec$Moderator]] / spec$Scale
    complete$moderator_scaled_centered <- scaled - mean(scaled)
    moderator_mean_original <- mean(complete[[spec$Moderator]])
    for (rho in c(0, 0.5, 0.8)) {
      fit <- fit_model(complete, rho)
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        Driver = spec$Driver, Nutrient = nutrient, Moderator = spec$Moderator,
        Moderator_unit = spec$Unit, rho = rho, effects = nrow(complete),
        studies = studies, experiments = experiments,
        moderator_mean = moderator_mean_original,
        intercept_at_mean = fit$beta[1], intercept_CI_low = fit$lo[1], intercept_CI_high = fit$hi[1],
        slope_lnRR = fit$beta[2], SE_model = fit$se_model[2],
        SE_cluster = fit$se_cluster[2], SE_used = fit$se_used[2],
        CI_low = fit$lo[2], CI_high = fit$hi[2], p_value = fit$p_value[2],
        percent_change_in_ratio_per_unit = 100 * expm1(fit$beta[2]),
        percent_CI_low = 100 * expm1(fit$lo[2]), percent_CI_high = 100 * expm1(fit$hi[2]),
        tau2_study = fit$tau[1], tau2_experiment = fit$tau[2], df = fit$df,
        convergence = fit$convergence,
        evidence_flag = if (studies < 5) "Very limited (<5 studies)" else "Exploratory"
      )
    }
    if (studies >= 4) {
      for (omitted in sort(unique(complete$Study_ID))) {
        reduced <- complete[complete$Study_ID != omitted, ]
        if (length(unique(reduced$Study_ID)) < 3) next
        fit <- fit_model(reduced, 0)
        loo_rows[[length(loo_rows) + 1L]] <- data.frame(
          Driver = spec$Driver, Nutrient = nutrient, Moderator = spec$Moderator,
          Moderator_unit = spec$Unit, omitted_study = omitted,
          studies = length(unique(reduced$Study_ID)), experiments = length(unique(reduced$Experiment_ID)),
          effects = nrow(reduced), slope_lnRR = fit$beta[2], CI_low = fit$lo[2],
          CI_high = fit$hi[2], p_value = fit$p_value[2],
          percent_change_in_ratio_per_unit = 100 * expm1(fit$beta[2]),
          convergence = fit$convergence
        )
      }
    }
  }
}

feasibility <- do.call(rbind, feasibility_rows)
results <- do.call(rbind, result_rows)
loo <- do.call(rbind, loo_rows)
write.csv(feasibility, file.path(out_dir, "environmental_contrast_meta_regression_feasibility.csv"), row.names = FALSE)
write.csv(results, file.path(out_dir, "environmental_contrast_meta_regression_results.csv"), row.names = FALSE)
write.csv(loo, file.path(out_dir, "environmental_contrast_meta_regression_leave_one_study_out.csv"), row.names = FALSE)

primary <- results[results$rho == 0, c("Driver", "Nutrient", "Moderator", "Moderator_unit",
                                      "studies", "experiments", "effects", "slope_lnRR",
                                      "CI_low", "CI_high", "p_value",
                                      "percent_change_in_ratio_per_unit", "percent_CI_low",
                                      "percent_CI_high", "evidence_flag")]
write.csv(primary, file.path(out_dir, "environmental_contrast_meta_regression_primary_summary.csv"), row.names = FALSE)

print(feasibility)
cat("\nPrimary rho=0 moderator estimates:\n")
print(primary)
