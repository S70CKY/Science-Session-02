# Libraries
library(readxl)
library(dplyr)

# Load and Prepare Data ----

df <- read_excel("data/ScienceSession04.xlsx", sheet = 2, skip = 2)

oldname <- colnames(df)
newname <- c("cohort", "Family ID", "child_id", "sample_id", "age_d", "age_m",
             "whz", "waz", "haz", "breast_milk", "formula", "solid", "diarrhea",
             "antibiotic_7d", "medication", "no_sequence", "serum_id", "barcode")
df2 <- df %>% rename_at(vars(oldname), ~newname)

df2 <- df2 %>% filter(!is.na(child_id))

data3 <- df2 %>%
  mutate(age_y = age_m / 12 + age_d / 365.25) %>%
  relocate(age_y, .after = age_m)

df_meta <- read_excel("data/ScienceSession04.xlsx", sheet = 1)
a <- df_meta %>% select("Family ID", Gender)
cleaned_data <- data3 %>% left_join(a, by = "Family ID", relationship = "many-to-many")

# Descriptive Statistics ----

# Summary
summary(cleaned_data)
summary(cleaned_data$whz)

fivenum(cleaned_data$whz)

# Minimum and maximum
min(cleaned_data$whz, na.rm = TRUE)
max(cleaned_data$whz, na.rm = TRUE)

# Mean and median
mean(cleaned_data$whz, na.rm = TRUE)
median(cleaned_data$whz, na.rm = TRUE)

# Standard deviation
sd(cleaned_data$whz, na.rm = TRUE)

# Interquartile range
quantile(cleaned_data$whz, probs = c(0.25, 0.75), na.rm = TRUE)
IQR(cleaned_data$whz, na.rm = TRUE)

# Histogram
hist(cleaned_data$age_y, main = "age_Y Score Distribution", xlab = "WHZ")

# Table of counts (cat)
table(cleaned_data$antibiotic_7d)

# Proportions (%)
prop.table(table(cleaned_data$antibiotic_7d))

# Crosstab between two categorical variables
table(cleaned_data$diarrhea, cleaned_data$antibiotic_7d)

# Inferential Statistics ----

# Q-Q plot
qqnorm(cleaned_data$whz, main = "Q-Q Plot of WHZ")
qqline(cleaned_data$whz, col = "red", lwd = 2)

# Shapiro-Wilk Test
shapiro.test(cleaned_data$waz[!is.na(cleaned_data$waz)])

# t-test:
t.test(whz ~ antibiotic_7d, data = cleaned_data)

# ANOVA
aov_result <- aov(haz ~ breast_milk + formula + solid, data = cleaned_data)
summary(aov_result)

# TukeyHSD

#Correlation test
cor.test(cleaned_data$age_y, cleaned_data$whz, use = "complete.obs")

cor.test(cleaned_data$age_y, cleaned_data$whz, method = "spearman", exact = FALSE)

# Chi-Square test
tab <- table(cleaned_data$diarrhea, cleaned_data$antibiotic_7d)
chisq.test(tab)

#Wilcoxon test
wilcox.test(whz ~ antibiotic_7d, data = cleaned_data)

#Kruskal-Wallis test
kruskal.test(haz ~ breast_milk, data = cleaned_data)
# dunn's test