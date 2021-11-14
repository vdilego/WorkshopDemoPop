
# --------------------------------------------------------------------------------------------- #
# Author: Vanessa di Lego and Markus Sauerberg
# --------------------------------------------------------------------------------------------- #

library(here)
library(data.table)
library(readxl)
library(tidyverse)
library(vroom)

# --------------------------------------------------------------------------------------------- #
# README
# --------------------------------------------------------------------------------------------- #
# this script already considers that the data has been approved for user from EU-SILC and that  #
# data has been manipulated using survey package, in order to retrieve health prevalence by age #
# --------------------------------------------------------------------------------------------- #

# if using Markus' original spreadsheet

hly.folder.silc<-here("HLY_SILC","Data")

dir.create(hly.folder.silc, showWarnings = FALSE, recursive = TRUE)

file.list <- list.files(here("WP5"),pattern='*.xlsx', recursive=T)

df <- list.files(path = here("WP5"),
                 full.names = T,
                 recursive = TRUE,
                 pattern = "*.xlsx") %>%
  tibble::as_tibble() %>%
  mutate(sheetName = map(value, readxl::excel_sheets)) %>%
  unnest(sheetName) %>%
  mutate(myFiles = purrr::map2(value, sheetName, function(x,y) {
    readxl::read_excel(x, sheet = paste(y))})) %>%
  unnest(myFiles) %>%
  select(!c(1,2))

saveRDS(df, file = file.path(hly.folder.silc, "HealthData.rds"))

# if using Markus' latest .txt files
# Read all the files that are .txt in the directory and create a FileName column to store filenames
file.list<- list.files(path = here("WP5"), recursive = TRUE,
                            pattern = "\\.txt$",
                            full.names = F)

df2 <- list.files(path = here("WP5"),
                 full.names = F,
                 recursive = TRUE,
                 pattern = "\\.txt$") %>%
  map(~ read.table(file.path(here("WP5"), .))) %>%
  reduce(rbind)

saveRDS(df2, file = file.path(hly.folder.silc, "HealthData_ext.rds"))
