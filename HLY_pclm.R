
# --------------------------------------------------------------------------------------------- #
# Author: Vanessa di Lego and Markus Sauerberg
# --------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------
#  Part 4: Healthy Life Years extrapolated by pclm
# -------------------------------------------------------------------------------------------

library(dplyr)
library(openxlsx)
library(eurostat)
library(HMDHFDplus)
# remotes::install_github("patrickaubert/healthexpectancies",ref='main')
library(healthexpectancies)
library(here)
library(data.table)
library(tidyverse)
library(styler)
library(forcats)
library(MortalityLaws)
library(purrr)
library(broom)
library(ggplot2)
library(gridExtra)
library(MortalitySmooth)
library(ungroup)
# loading useful functions into environment
source(here("0_Functions.R"))

# setting path and loading health data
hly.folder.silc <- here("HLY_SILC", "Data")

hly.data_fra <- readRDS(file.path(hly.folder.silc, "HealthData_ext.rds")) %>%
  mutate(Country=as.factor(Country) %>%
           fct_recode("GB"="UK")) %>%
  select(1:7) %>%
  relocate(Country, .before = Year) %>%
  mutate(Sex=as.factor(Sex) %>%
           fct_recode("Female"="F", "Male"="M")) %>%
  filter(Country=="FR" & Year==2018 & Sex=="Female" & Age<=81) %>%
  rename(Limited=Limited.weighted,
         Healthy=Unlimited.weighted,
         Prev=Prev.weighted) %>%
  mutate(Prev=coalesce(Prev, ((sum(Limited[Age<=19], na.rm = T)/
                                (sum(Limited[Age<=19], na.rm = T)+sum(Healthy[Age<=19],
                                                                      na.rm = T))))/2))



  hly.data_fra[is.na(hly.data_fra)] <- 0

options(scipen=999)

# install.packages("devtools")

library(devtools)
# requires the development version of rstan, sorry!
install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
install_github("timriffe/DemoTools", dependencies=T)
library(DemoTools)


## adding new column to data.table with extrapolated rates

# use this for all countries  graduation.
hly.pclm <- function(x,y){

  # last bin interval
  nlast  <- 19

  #  exposure - here it will be the total of limited+unlimited
  #offset <- z
  Prev<-y

  M2      <- pclm(x = x,
                  y = Prev,
                  nlast = nlast,
                  #offset = z,
                  control = list(lambda = 1))$fitted
  Age<-0:99
  df <- data.frame(Age,M2)
  colnames(df) <- c("Age","Prev")
  return(df)

}


#setting data.table
setDT( hly.data_fra)
hly.data_fra$Age<-as.numeric(hly.data_fra$Age)
hly.data_fra$Prevalence<-as.numeric(hly.data_fra$Prev)

# takin out hideous scientific notation at reporting rates
options(scipen=999)

x<-seq(0,81,by=1)
#apply function for estimating pclm into single ages
prev.grad<-hly.data_fra[, hly.pclm(x,Prev)]
setDT(prev.grad)
# adding new column to data.table with extrapolated rates

prev.grad<-prev.grad[ ,Prev_extrap:= lt_rule_m_extrapolate(Prev,Age,x_fit = 60:75,
                                                        x_extr = 75:99, law="kannisto_makeham" )$values]

ggplot(prev.grad, aes(Age,Prev_extrap))+
  geom_line()

