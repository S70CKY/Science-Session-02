# Installation ----
install.packages("")

# Library ----
library(tidyverse)
library(readr)
library(readxl)

# Basics ----
## Data Types ----
a <- 3.14
b <- "Hello"
c <- TRUE

d <- factor(c("low", "medium", "high"))

## Data Structure ----
vec <- c(1, 2, 3)
m <- matrix(c(10, 20, 30, 40), nrow = 2, ncol = 2)
lst <- list(name = "Alice", age = 30, passed = TRUE)
df <- data.frame(id = 1:3, score = c(88, 95, 70))
tb <- tibble(id = 1:3, score = c(88, 95, 70))

# Import ----
test_text1 <- read_delim("data/test_txt1.txt")
test_csv <- read_csv("data/test_csv.csv")
test_xls <- read_xls("data/test_xls.xls")
test_xlsx <- read_xlsx("data/test_xlsx.xlsx")
test_text1 <- read_tsv("data/test_txt1.txt")

# data1 <- read_excel("data/health_data_1.xlsx")
# data2 <- read_excel("data/health_data_2.xlsx")
# data3 <- read_excel("data/health_data_3.xlsx")

# Tidying ----
untidy_data <- tibble(
  name = c("Ana","Bob","Cara"),
  meds = c("advil 600mg 2xday","tylenol 650mg 4xday", "advil 200mg 3xday")
)
untidy_data

## Transpose ----
scores <- tibble(
  student = c("John", "Jane"),
  math = c(80, 90),
  science = c(85, 88)
)

long_scores <- scores %>%
  pivot_longer(cols = c(math, science), names_to = "subject", values_to = "score")

wide_again <- long_scores %>%
  pivot_wider(names_from = subject, values_from = score)

## Demonstration ----
demo_data <- read_csv("data/yrbss_demo.csv")

glimpse(demo_data)
str(demo_data)       
head(demo_data)
summary(demo_data)
class(demo_data)

data <- demo_data %>%
  mutate(sex = case_when(
    sex == "Male" ~ 1,
    sex == "Female" ~ 0,
    TRUE ~ NA_real_ 
  ))

### rename ----
demo_data %>% rename(id = record)

### select ----
demo_data %>% filter(bmi < 5) 

### filter ----
demo_data %>% filter(race4 == "Asian")

### mutate ----
newdata <- demo_data %>% 
  mutate(height_m = sqrt(stweight / bmi))

demo_data2 <- demo_data %>%
  mutate(
    bmi_group = case_when(
      bmi < 18.5 ~ "underweight",             
      bmi >= 18.5 & bmi <= 24.9 ~ "normal",
      bmi > 24.9 & bmi <= 29.9 ~ "overweight",
      bmi > 29.9 ~ "obese")
    )

### separate ----
new <- demo_data %>% 
  separate(age,c("1","2","3","4","5"),
           sep = " ") %>%
  select("1":"5")

### unite ----
new2 <- demo_data %>% 
  unite("sexgr", sex, grade, sep=":") %>%
  select(sexgr)

### remove rows with missing data ----
demo_data %>% na.omit()

### arrange ----
a <- demo_data %>% arrange(desc(bmi))

### summarize ----
glimpse(demo_data)
summary(demo_data)
summary(newdata)

# Join (merge) ----
df1 <- tibble(
  id = c(1, 2, 3),
  name = c("Alice", "Bob", "Charlie")
)

df2 <- tibble(
  id = c(2, 3, 4),
  score = c(88, 92, 75)
)
df1
df2

left_join(df1, df2, by = "id")
inner_join(df1, df2, by = "id")
full_join(df1, df2, by = "id")
right_join(df1, df2, by = "id")
anti_join(df1, df2, by = "id")

# höo

# Export ----
write_csv(demo_data, "data/newdata.csv")

untidy_data %>% 
  separate(col = meds, into = c("med_name","dose_mg","times_per_day"), sep=" ") %>% 
  mutate(times_per_day = as.numeric(str_remove(times_per_day, "xday")),         
         dose_mg = as.numeric(str_remove(dose_mg, "mg")))
