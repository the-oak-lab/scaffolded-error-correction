library(tidyverse)
library(lme4)
library(ggeffects)
library(broom)
library(patchwork)

d <- read_csv("anonymized_data.csv")

# Contrasts ---------------------------------------------------------------

d$explanation_v_cr <- ifelse(d$condition == "correct_answer", -.75, .25)

d$socratic_vs_othexp <- case_when(
  d$condition == "socratic" ~ .66667,
  grepl("explanation", d$condition) ~ -.3333,
  TRUE ~ 0
)

d$personalized_vs_genexp <- case_when(
  d$condition == "personalized_explanation" ~ .5,
  d$condition == "generic_explanation" ~ -.5,
  TRUE ~ 0
)


# Standardizing pretest scores --------------------------------------------

d$pct_pretest_z <- d$pct_pretest %>% scale() %>% as.numeric()


# Labels and colors for graphs --------------------------------------------

condition_colors <- c("#B0B4B8", "#8FB4D4", "#2C5F8A", "#E8A33D")

d$condition_str <- factor(
  d$condition,
  c("correct_answer", "generic_explanation", "personalized_explanation", "socratic"),
  c("Correct\nAnswer", "Generic\nExplanation", "Personalized\nExplanation", "Socratic\nTutoring")
)

d$condition_str_linegraph <- factor(
  d$condition,
  c("correct_answer", "generic_explanation", "personalized_explanation", "socratic"),
  c("Correct Ans.", "Generic", "Personalized", "Socratic")
)


# Descriptive Stats -------------------------------------------------------


corAPA <- function(x, missingMethod=c("pairwise","complete"),
                   format = c("kable", "df"),
                   descriptives = T,
                   removeTriangle=c("upper", "lower"),
                   decimals=2, cutoffs=c(.001, .01, .05)){
  
  df <- x
  format = format[1]
  #Compute correlation matrix
  x <- as.matrix(x)
  correlation_matrix<-psych::corr.test(x, use = missingMethod[1])
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$p # Matrix of p-value
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < cutoffs[1], "***", ifelse(p < cutoffs[2], "** ", ifelse(p < cutoffs[3], "*  ", ifelse(p < .1, "    ", "    "))))
  
  ## Define notions for significance levels; spacing is important.
  mystars <- ifelse(p < cutoffs[1], "***", ifelse(p < cutoffs[2], "** ", ifelse(p < .05, "*  ", ifelse(p < .1, "    ", "    "))))
  
  ## trunctuate the correlation matrix to number of decimals
  R <- format(round(cbind(rep(-1.11, ncol(x)), R), decimals))[,-1]
  
  ## build a new matrix that includes the correlations with their apropriate stars
  Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
  diag(Rnew) <- paste(diag(R), " ", sep="")
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep="")
  
  ## remove upper triangle of correlation matrix
  if(removeTriangle[1]=="upper"){
    Rnew <- as.matrix(Rnew)
    Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew) %>% tibble::rownames_to_column(var = "measure")
  }
  
  ## remove lower triangle of correlation matrix
  if (removeTriangle[1]=="lower"){
    Rnew <- as.matrix(Rnew)
    Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
    Rnew <- as.data.frame(Rnew) %>% tibble::rownames_to_column(var = "measure")
  }
  
  ## add numbering
  Rnew <- Rnew %>% dplyr::mutate(measure = stringr::str_c(dplyr::row_number(), ". ", measure))
  
  ## add descriptives
  if(descriptives){
    msd <- df %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = "name") %>%
      dplyr::group_by(name) %>%
      dplyr::summarise(mean = mean(value, na.rm = T),
                       sd = sd(value, na.rm = T)) %>%
      tidyr::pivot_longer(-name, names_to = "measure") %>%
      tidyr::pivot_wider(names_from = name, values_from = value) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, decimals))) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), as.character))
    
    Rnew <- dplyr::bind_rows(Rnew, msd)
  }
  
  ## return the correlation matrix
  if(format == "kable"){
    knitr::kable(Rnew, format = 'html', booktabs = T) %>% kableExtra::kable_styling()
  }
  
  else if(format == "df"){
    return(Rnew)
  }
  
}


d %>% select(pct_pretest, n_practice_correct, pct_practice_aligned:frustrating_f) %>%
  corAPA()

d %>%
  group_by(condition) %>%
  summarise(across(c(pct_pretest, n_practice_correct:frustrating_f),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}")) %>%
  pivot_longer(-condition) %>%
  mutate(score = ifelse(grepl("mean", name), "mean", "sd")) %>%
  mutate(name = str_remove(name, "_mean") %>% str_remove("_sd")) %>%
  pivot_wider(names_from = score, values_from = value) %>%
  pivot_wider(names_from = condition, values_from = mean:sd, names_vary = "slowest") %>%
  knitr::kable() %>% kableExtra::kable_styling()

# Primary Models ----------------------------------------------------------

primary_outcomes <- c(
  "pct_practice_aligned",
  "pct_pretest_aligned",
  "pct_challenge",
  "si_session_f",
  "int_topic_f",
  "practice_minutes",
  "jol_f",
  "feedback_eff_f",
  "cost_f",
  "frustrating_f"
)

# Preregistered model
primary_rhs <- "(explanation_v_cr + socratic_vs_othexp + personalized_vs_genexp) * pct_pretest_z"

# Fit all models

primary_models      <- list()
primary_coef_tables <- list()
primary_fit_tables  <- list()

for (y in primary_outcomes) {
  
  m <- lm(as.formula(paste(y, "~", primary_rhs)), data = d)
  
  # Scaling constants: grand mean, and residual SD after removing condition
  mean_y <- mean(d[[y]], na.rm = TRUE)
  sd_y   <- sigma(lm(as.formula(paste(y, "~ condition")), data = d))
  
  # Coefficients, raw and standardized
  coefs <- tidy(m, conf.int = TRUE) %>% select(term:p.value)
  
  shift <- ifelse(coefs$term == "(Intercept)", mean_y, 0)   # only the intercept is re-centered
  
  coefs$std_estimate <- (coefs$estimate - shift) / sd_y
  coefs$outcome      <- y
  
  # Model-level fit stats
  fit <- glance(m)
  fit$n       <- nobs(m)
  fit$outcome <- y
  
  # Save everything under the name of the outcome
  primary_models[[y]]      <- m
  primary_coef_tables[[y]] <- coefs
  primary_fit_tables[[y]]  <- fit
}

rm(y, m, mean_y, sd_y, coefs, shift, fit)   

# Combine and print

primary_summaries <- bind_rows(primary_coef_tables) %>%
  select(outcome, term, estimate, std_estimate, std.error, statistic, p.value)

primary_fit <- bind_rows(primary_fit_tables) %>%
  select(outcome, n, everything())

primary_summaries %>%
  knitr::kable(digits = 3, format = "html") %>%
  kableExtra::kable_styling(full_width = FALSE)


# Plots -------------------------------------------------------------------

measure_labels <- c(
  pct_practice_aligned = "Practice-Aligned Questions",
  pct_pretest_aligned  = "Pretest-Aligned Question",
  pct_challenge        = "Challenge Question"
)


# Performance Outcomes ----------------------------------------------------

d %>%
  pivot_longer(all_of(names(measure_labels)),
               names_to = "measure", values_to = "pct") %>%
  mutate(measure = factor(measure,
                          levels = names(measure_labels),
                          labels = measure_labels)) %>%
  ggplot(aes(x = condition_str, y = pct, fill = condition_str)) +
  stat_summary(geom = "bar", fun = mean) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = .25) +
  facet_wrap(~ measure, nrow = 1) +
  coord_cartesian(ylim = c(0, .8)) +
  scale_y_continuous(breaks = seq(0, 1, by = .2), expand = c(0, 0),
                     labels = scales::percent) +
  scale_x_discrete(labels = scales::label_wrap(10)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Percent Correct") +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(size = rel(0.7))
  )

ggsave("fig3_performance.pdf",
       width = 190, height = 80, units = "mm")


# Cost and frustration ----------------------------------------------------

p_cf_a <- d %>%
  ggplot(aes(x = condition_str, y = cost_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 6)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Perceived Cost") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )

p_cf_b <- d %>%
  ggplot(aes(x = condition_str, y = frustrating_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 6)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Frustration") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )


# Cost interaction --------------------------------------------------------

m_cost_graph <- lm(cost_f ~ condition_str_linegraph * pct_pretest, data = d)

pred_cost <- ggpredict(
  m_cost_graph,
  terms = c("pct_pretest [0:0.67 by=0.01]", "condition_str_linegraph")
)

p_cf_c <- ggplot(pred_cost, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  coord_cartesian(ylim = c(1, 6)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 1, by = .2), expand = c(0, 0),
                     labels = scales::percent) +
  scale_color_manual(values = condition_colors) +
  scale_fill_manual(values = condition_colors) +
  labs(x = "Pretest score", y = "Perceived Cost",
       color = NULL, fill = NULL) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.text = element_text(size = rel(0.7)),
    legend.background = element_rect(fill = alpha("white", 0.8),
                                     color = "black", linewidth = .2)
  )

p_cf_a + p_cf_b + p_cf_c + plot_layout(ncol = 3)

ggsave("fig4_cost_frustration.pdf",
       width = 190, height = 80, units = "mm")


# Interest measures -------------------------------------------------------

p_int_1 <- d %>%
  ggplot(aes(x = condition_str, y = si_session_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 6)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Situational Interest") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )

p_int_2 <- d %>%
  ggplot(aes(x = condition_str, y = int_topic_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 5)) +
  scale_y_continuous(breaks = seq(1, 5, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Interest in Constraint Equations") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )

p_int_1 + p_int_2 + plot_layout(ncol = 2)

ggsave("fig5_interest.pdf",
       width = 130, height = 80, units = "mm")


# JOL and Effectiveness ---------------------------------------------------

p_jol_eff_1 <- d %>%
  ggplot(aes(x = condition_str, y = jol_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 5)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Judgment of Learning") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )

p_jol_eff_2 <- d %>%
  ggplot(aes(x = condition_str, y = feedback_eff_f, fill = condition_str)) +
  geom_bar(stat = "summary") +
  geom_errorbar(stat = "summary", width = .25) +
  coord_cartesian(ylim = c(1, 6)) +
  scale_y_continuous(breaks = seq(1, 6, by = 1), expand = c(0, 0)) +
  scale_fill_manual(values = condition_colors) +
  labs(x = NULL, y = "Perceived Effectiveness of Feedback") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = rel(0.7))
  )

p_jol_eff_1 + p_jol_eff_2 + plot_layout(ncol = 2)

ggsave("fig6_jol_effectiveness.pdf",
       width = 130, height = 80, units = "mm")


# Process Models in the Socratic Condition --------------------------------

d$n_practice_correct_c <- d$n_practice_correct - mean(d$n_practice_correct)

process_outcomes <- c(
  "pct_practice_aligned",
  "jol_f",
  "si_session_f",
  "int_topic_f",
  "cost_f",
  "frustrating_f"
)

process_rhs <- "n_practice_fixed + n_practice_correct_c"

# Standardized outcomes: grand mean, and residual SD after removing condition
for (y in process_outcomes) {
  sd_y <- sigma(lm(reformulate("condition", response = y), data = d))
  d[[paste0(y, "_z")]] <- (d[[y]] - mean(d[[y]], na.rm = TRUE)) / sd_y
}

rm(y, sd_y)   # clear loop scratch variables

soc <- d %>% filter(condition == "socratic")

process_models      <- list()
process_coef_tables <- list()
process_test_tables <- list()

for (y in process_outcomes) {
  
  y_z <- paste0(y, "_z")
  
  m <- lm(as.formula(paste(y_z, "~", process_rhs)), data = soc)
  
  # Coefficients (outcome is already standardized)
  coefs <- tidy(m, conf.int = TRUE) %>%
    mutate(var = y)
  
  # Linear hypothesis: are the two slopes equal?
  lh <- car::linearHypothesis(m, "n_practice_fixed = n_practice_correct_c")
  
  tests <- tibble(
    var    = y,
    F      = lh$F[2],
    df1    = lh$Df[2],
    df2    = lh$Res.Df[2],
    p_diff = lh$`Pr(>F)`[2]
  )
  
  # Save under this outcome's name
  process_models[[y]]      <- m
  process_coef_tables[[y]] <- coefs
  process_test_tables[[y]] <- tests
}

rm(y, y_z, m, coefs, lh, tests) 

process_coefs <- bind_rows(process_coef_tables)
process_tests <- bind_rows(process_test_tables)

process_coefs %>%
  knitr::kable(digits = 3, format = "html") %>%
  kableExtra::kable_styling(full_width = FALSE)

# Plot

process_coefs$reverse <- process_coefs$var %in% c("frustrating_f", "cost_f")

process_plot_df <- process_coefs %>%
  filter(term != "(Intercept)") %>%
  mutate(
    est = if_else(reverse, -estimate,  estimate),
    lo  = if_else(reverse, -conf.high, conf.low),
    hi  = if_else(reverse, -conf.low,  conf.high),
    term = factor(
      term,
      levels = c("n_practice_fixed", "n_practice_correct_c"),
      labels = c("Benefit of Resolving an Error",
                 "Benefit of a Correct First Attempt")
    ),
    label = factor(
      var,
      c("jol_f", "int_topic_f", "si_session_f",
        "frustrating_f", "cost_f", "pct_practice_aligned"),
      c("Judgment of learning", "Interest in constraint equations",
        "Situational interest", "Frustration (reversed)",
        "Perceived cost (reversed)", "Near-transfer performance")
    )
  )

ggplot(process_plot_df, aes(x = est, y = label, color = term, shape = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi),
                  position = position_dodge(width = 0.6), size = 0.5) +
  scale_color_manual(values = c("#E8A33D", "#2C5F8A")) +
  labs(x = "Standardized Units", y = NULL, color = NULL, shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  guides(color = guide_legend(reverse = TRUE),
         shape = guide_legend(reverse = TRUE))

ggsave("fig7_resolved_errors.pdf", width = 190, height = 100, units = "mm")
