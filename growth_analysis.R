# growth_analysis.R
# Tree circumference growth analysis from growth_data.csv

library(dplyr)
library(tidyr)
library(ggplot2)

# ---- Load data ----
growth_data <- read.csv("growth_data.csv", stringsAsFactors = FALSE)

# 6. Column names
print(colnames(growth_data))

# ---- Reshape data ----
growth_long <- growth_data %>%
  pivot_longer(
    cols = starts_with("Circumf"),
    names_to = "Year",
    values_to = "Circumference"
  ) %>%
  mutate(Year = as.numeric(gsub("Circumf_|_cm", "", Year)))

# 7. Mean and SD at start/end
start_year <- min(growth_long$Year, na.rm = TRUE)
end_year   <- max(growth_long$Year, na.rm = TRUE)

growth_summary <- growth_long %>%
  filter(Year %in% c(start_year, end_year)) %>%
  group_by(Site, Year) %>%
  summarise(
    mean_circ = mean(Circumference, na.rm = TRUE),
    sd_circ   = sd(Circumference, na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  )
print(growth_summary)

# 8. Box plot start vs end
plot_data <- growth_long %>%
  filter(Year %in% c(start_year, end_year)) %>%
  mutate(Timepoint = factor(Year, levels = c(start_year, end_year),
                            labels = c(paste0("Start_", start_year), paste0("End_", end_year))))

p <- ggplot(plot_data, aes(x = Timepoint, y = Circumference, fill = Site)) +
  geom_boxplot(position = position_dodge(0.75)) +
  theme_minimal() +
  labs(title = paste0("Tree Circumference: Start (", start_year, ") vs End (", end_year, ")"),
       x = "Timepoint", y = "Circumference (cm)")
print(p)

ggsave("circumference_start_end_boxplot.png", p, width = 8, height = 5)

# 9. Mean growth over last 10 years (per-tree)
year_before_10 <- end_year - 10

end_df    <- growth_long %>% filter(Year == end_year)
before_df <- growth_long %>% filter(Year == year_before_10)

paired <- inner_join(end_df, before_df, by = c("TreeID", "Site"), suffix = c("_end", "_before"))
paired <- paired %>%
  mutate(growth_10yr = Circumference_end - Circumference_before)

mean_growth_by_site <- paired %>%
  group_by(Site) %>%
  summarise(
    mean_growth_10yr = mean(growth_10yr, na.rm = TRUE),
    sd_growth_10yr   = sd(growth_10yr, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(mean_growth_by_site)

write.csv(mean_growth_by_site, "mean_10yr_growth_by_site.csv", row.names = FALSE)

# 10. t-test
ttest_result <- t.test(growth_10yr ~ Site, data = paired)
print(ttest_result)
capture.output(ttest_result, file = "ttest_10yr_growth_result.txt")
