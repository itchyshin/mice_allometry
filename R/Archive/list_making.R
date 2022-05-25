# preparing dat_list2

library(tidyverse)
library(here)

allometryNew <- readRDS(here("data", "allometryNEW.rds"))

dat_list2 <- split(allometryNew, allometryNew$parameter_name)

saveRDS(dat_list2, file= "dat_list2.rds")