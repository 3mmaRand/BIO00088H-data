# load iris dataset
data(iris)
# what are the columns
names(iris)
# what is the range of values for each column
summary(iris)
# what is the mean of each column
colMeans(iris)
# the mean using the tidyverse
library(tidyverse)
iris %>% 
  summarise(across(everything(), mean))
# the mean for each species
iris %>% 
  group_by(Species) %>% 
  summarise(across(everything(), mean))
# the mean, standard deviation and length for each species
iris %>% 
  group_by(Species) %>% 
  summarise(across(everything(), list(mean, sd, length)))

# write a function to calculate the standard error
# of the mean
se <- function(x) {
  sd(x) / sqrt(length(x))
}








