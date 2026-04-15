################################################################################
# Scatter plot: Consensus (X_bar) vs Disagreement (S) for each ordinal feature
# Output: fig_scatter_cd.pdf (for inclusion in manuscript)
################################################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(forcats)
  library(readxl)
})

# ── Data Loading ──────────────────────────────────────────────
data_path <- "../../data/data.xlsx"
data <- read_excel(data_path)

data_clean <- data %>%
  slice(-(1:2)) %>%
  rename(Pre_DWI_ADC = `DWI-ADC...12`, Post_DWI_ADC = `DWI-ADC...23`, pCR = `...28`) %>%
  mutate(across(c(Pre_EMVI, normalized_wall, ulcer, split_scar,
                  Pre_cT, Pre_DWI, Pre_ADC, Post_cT, Fibrosis_grading,
                  Post_DWI, Post_ADC, Post_restriction, MERCURY, ESGAR,
                  pCR, Center), as.numeric),
         reviewer = as.numeric(reviewer))

ordinal_vars <- c("Pre_cT", "Pre_DWI", "Pre_ADC", "Post_cT", "Fibrosis_grading",
                  "Post_DWI", "Post_ADC", "Post_restriction", "MERCURY", "ESGAR")

# Labels matching paper (Table 6 order by |r|)
var_labels <- c(
  "Pre_DWI" = "Pre_DWI",
  "Post_DWI" = "Post_DWI",
  "Pre_ADC" = "Pre_ADC",
  "ESGAR" = "ESGAR",
  "Pre_cT" = "Pre_cT",
  "Fibrosis_grading" = "Fibrosis",
  "Post_ADC" = "Post_ADC",
  "Post_restriction" = "Post_restr.",
  "Post_cT" = "Post_cT",
  "MERCURY" = "MERCURY"
)

cd_data <- data_clean %>%
  group_by(SubjectNo) %>%
  summarise(
    pCR = first(pCR), Center = first(Center),
    across(all_of(ordinal_vars), list(
      cons = ~ mean(.x, na.rm = TRUE),
      disp = ~ sd(.x, na.rm = TRUE)
    )),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_disp"), ~ replace_na(.x, 0)))

# ── Build long-format data for faceted plot ───────────────────
plot_data <- map_dfr(ordinal_vars, function(v) {
  x <- cd_data[[paste0(v, "_cons")]]
  s <- cd_data[[paste0(v, "_disp")]]
  r2_lin  <- summary(lm(s ~ x))$r.squared
  r2_quad <- summary(lm(s ~ x + I(x^2)))$r.squared
  label <- paste0(var_labels[v],
                  "\n(lin R\u00b2=", sprintf("%.2f", r2_lin),
                  ", quad R\u00b2=", sprintf("%.2f", r2_quad), ")")
  tibble(
    Consensus = x,
    Disagreement = s,
    Feature = label,
    r2_quad = r2_quad
  )
})

# Order by quadratic R² descending
plot_data <- plot_data %>%
  mutate(Feature = fct_reorder(Feature, -r2_quad))

# ── Create scatter plot ───────────────────────────────────────
p <- ggplot(plot_data, aes(x = Consensus, y = Disagreement)) +
  geom_point(alpha = 0.4, size = 1.2, color = "steelblue") +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = TRUE, color = "red3", linewidth = 0.7, alpha = 0.2) +
  facet_wrap(~ Feature, scales = "free", ncol = 5) +
  labs(
    x = expression(bar(X)[k] ~ "(consensus: mean across readers)"),
    y = expression(S[k] ~ "(disagreement: SD across readers)")
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text = element_text(size = 7, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

# ── Save ──────────────────────────────────────────────────────
out_path <- "../../figures/fig_scatter_cd.pdf"
ggsave(out_path, p, width = 10, height = 4.5, dpi = 300)
cat("Saved:", out_path, "\n")
