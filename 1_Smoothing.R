# --------------------------------------------------------------------------------------------- #
# Author: Vanessa di Lego and Markus Sauerberg
# --------------------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------------------
#  Part 1: Healthy Life Years and different smoothing methods for health prevalence by age
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
#devtools::install_github("karthik/wesanderson")
library(wesanderson)
library(ggthemes)
# devtools::install_github("alburezg/suffrager")
library(suffrager)
library(remotes)
#remotes::install_github("cran/MortalitySmooth", dependencies = T)
library(MortalitySmooth)

# loading useful functions into environment
source(here("0_Functions.R"))
# setting path and loading health data
hly.folder.silc <- here("HLY_SILC", "Data")

hly.data <- readRDS(file.path(hly.folder.silc, "HealthData_ext.rds")) %>%
  mutate(Country=as.factor(Country) %>%
           fct_recode("GB"="UK")) %>%
  select(1:7) %>%
  relocate(Country, .before = Year) %>%
  mutate(Sex=as.factor(Sex) %>%
           fct_recode("Female"="F", "Male"="M"))


### Download HMD data using MortalityLaws Package. Add your username and password
# select countries of interest. We will select only those that we can match with the health data from
# EU-SILC and that are comparable.

countries<-c("DEUTNP","ESP","FRATNP",
             "ITA","GBR_NP","POL","DNK")

deaths <- ReadHMD(what="Dx",
                  countries=countries,
                  interval="1x1",
                  username = "username here",
                  password = "your password here",
                  save=F)$data %>%
  pivot_longer(4:6, names_to="sex", values_to="deaths")

exposure <- ReadHMD(what="Ex",
                    countries=countries,
                    interval="1x1",
                    username = "username here",
                    password = "your password here",
                    save=F)$data %>%
  pivot_longer(4:6, names_to="sex", values_to="exposure")

# Select year you wish to analyse. We will concentrate on years after 2014 due to health prevalence
# from EU-SILC and allow cross-year comparability.

all<-left_join(deaths,exposure,by=c("country","Year","Age","sex")) %>%
  filter(sex%in%c("Female","Male")
         & Year%in%2004:2019) %>%
  arrange(country,Year,sex,Age)

# First, we close the HMD life tables at age 80+. I kept Germany but estimates are not valid for HLY yet, since
# they stop at age 75. In the next analysis where I extrapolate health at older ages Germany is ok.
# In this case, I just take Germany out of the graphs later on.

# To close the HMD life tables at age 80+ we sum deaths and exposures above age 80, compute the truncated death rates
# at age 80 (nmx), and apply Markus' lifetable function.

lt_80<-all %>%
  group_by(country,Year,sex) %>%
  mutate(deaths=case_when(Age>=80~sum(deaths[Age>=80]),TRUE~deaths),
         exposure=case_when(Age>=80~sum(exposure[Age>=80]),TRUE~exposure),
         nmx=deaths/exposure) %>%
  filter(Age<=80) %>%
  ungroup() %>%
  arrange(country,Year,sex,Age) %>%
  group_by(country,Year,sex) %>%
  group_modify(~life.table(.x$nmx), .keep=T) %>%
  mutate(Age=0:80) %>%
  relocate(Age, .after = Year) %>%
  select(1:4,8,10,12) %>%
  rename(Country=country,Sex=sex)

# renaming the levels and labels of country factor to match the health database. I am also
# switching everything to ISO2 codes in case we make a map some day or want to use country.code package,
# but this is not necessary.

lt_80$Country <- factor(lt_80$Country,
                        labels=c("DE","DK","ES","FR",
                                 "IT","PL","GB"),
                        levels=c("DEUTNP","DNK","ESP",
                                 "FRATNP","ITA","POL",
                                 "GBR_NP"))
#saveRDS(lt_80, file = file.path(hly.folder.silc, "lifetable_80+.rds"))

# Join health data with mortality and group by age - Note: health data goes only until 80+ (in case of Germany to 75+).
# First we perform the analysis considering 80+ and then we employ other methods to extend the health prevalence
# for older ages. Groupping age only by 5 years as it´s more used.


hly.80<-hly.data %>%
  filter(Age<81) %>%
  left_join(lt_80,
            by=c("Country","Year","Age","Sex")) %>%
  rename(Limited=Limited.weighted,
         Healthy=Unlimited.weighted,
         Prevalence=Prev.weighted) %>%
  group_by(Country,Year,Sex) %>%
  mutate(AgeCat5=cut(Age, c(0,1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,Inf),
                     c("0-1","1-5","5-10","10-15","15-20",
                       "20-25","25-30", "30-35","35-40",
                       "40-45","45-50","50-55","55-60",
                       "60-65","65-70","70-75","75-80","80+"), include.lowest=TRUE, right=F),
         Lx.5=case_when(AgeCat5=="0-1"~ Lx[Age=="0"],TRUE~ 0)) %>%
  ungroup() %>%
  group_by(Country,Year,Sex,AgeCat5) %>%
  mutate(Lx.5=case_when(Lx.5==0~ sum(Lx),TRUE~ Lx.5),
         Limited.5=sum(Limited),
         Healthy.5=sum(Healthy),
         Prev.5=(Limited.5/(Limited.5+Healthy.5)),
         Prev.5=case_when(AgeCat5=="15-20"~ (sum(Limited[Age<=19], na.rm = T)/
                                               (sum(Limited[Age<=19], na.rm = T)+sum(Healthy[Age<=19],
                                                                                     na.rm = T))),TRUE~Prev.5))%>%
  ungroup() %>%
  group_by(Country,Year,Sex) %>%
  mutate(Prev.5=coalesce(Prev.5,(unique(Prev.5[AgeCat5=="15-20"])/2)),
         Prev.1=coalesce(Prevalence,0),
         Prev.1_half=case_when(Prev.1==0~ Prev.5,TRUE~ Prev.1))


#hly.no_disability0[is.na(hly.no_disability0)] <- 0

# estimating smoothed prevalences by single age and 5 year age group, using splines and polynomial fit

Health.table.80 <- hly.80 %>%
  group_by(Country,Year,Sex) %>%
  # no disability before age 17
  mutate(Prev.1_spline_0=case_when(Prev.1>0~
                                     smooth.spline(Prev.1)$y, TRUE ~ 0),
         #prevalence before age of 15 years is half of the prevalence of the next age interval (EU assumption).
         Prev.1_spline_half= smooth.spline(Prev.1_half)$y,
         Prev.1_penalty=smooth.spline(Prev.1_half, spar=0.7)$y,
         Prev.1_poly=prevalence_to_polynomial(Prev.1_half, agemin = 0,agemax = 80),
         # by 5-year age interval
         #prevalence before age of 15 years is half of the prevalence of the next age interval (EU assumption).
         Prev.5_spline_half=smooth.spline(Prev.5)$y,
         Prev.5_penalty=smooth.spline(Prev.5, spar=0.7)$y,
         Prev.5_poly=prevalence_to_polynomial(Prev.5, agemin = 0,agemax = 80))

View(Health.table.80)

# HLY for smoothed and unsmoothed prevalence, ending at 80+ single ages first

HLY.1 <- Health.table.80 %>%
  group_by(Country,Year,Sex) %>%
  # single ages
  mutate(LE=rev(cumsum(rev(Lx)))/lx,

         Lx.healthy.0=Lx*(1-Prev.1),
         Lx.healthy_half=Lx*(1-Prev.1_half),
         Lx.healthy.spline_0=Lx*(1-Prev.1_spline_0),
         Lx.healthy.spline_half=Lx*(1-Prev.1_spline_half),
         Lx.healthy.penalty=Lx*(1-Prev.1_penalty),
         Lx.healthy.poly=Lx*(1-Prev.1_poly),

         HLY.1=rev(cumsum(rev(Lx.healthy.0)))/lx,
         HLY.1_half=rev(cumsum(rev(Lx.healthy_half)))/lx,
         HLY.1_spline.0=rev(cumsum(rev(Lx.healthy.spline_0)))/lx,
         HLY.1_spline.half=rev(cumsum(rev(Lx.healthy.spline_half)))/lx,
         HLY.1_penalty=rev(cumsum(rev(Lx.healthy.penalty)))/lx,
         HLY.1_poly=rev(cumsum(rev(Lx.healthy.poly)))/lx)

# create subset to check differences

HLY_select<-HLY.1 %>%
  filter(Country%in% c("FR","DK") & Age%in%c("0","65")&Year==2017) %>%
  mutate(Prop.healthy=(HLY.1_half/LE)*100,
         Prop.healthy_spline=(HLY.1_spline.half/LE)*100,
         Prop.healthy_penalty=(HLY.1_penalty/LE)*100,
         Prop.healthy_poly=(HLY.1_poly/LE)*100)

write.table(HLY_select, "hly.csv",sep=",",row.names = F)

# create long dataset to perform graphs for single age ending in 80+

prev.compare<-HLY.1 %>%
  group_by(Country,Year,Sex,Age) %>%
  pivot_longer(starts_with("Prev.1"), names_to="Type",values_to="Value",
               names_repair = "unique")


prev.compare$Country<-as.factor(prev.compare$Country)
prev.compare$Country<- factor(prev.compare$Country,
                              levels=c("DE", "DK","ES","FR", "IT","PL","GB"),
                              labels = c("Germany","Denmark","Spain","France", "Italy","Poland","UK"))

prev.compare$Type<-as.factor(prev.compare$Type)
prev.compare$Type<- factor(prev.compare$Type,
                           levels=c("Prev.1","Prev.1_half","Prev.1_penalty",
                                    "Prev.1_poly","Prev.1_spline_0","Prev.1_spline_half"),
                           labels = c("No smoothing","No smoothing\n half disability",
                                      "Penalty","Polynomial", "Spline\nno Disability","Spline half\nthe disability"))


X11()
ggplot(prev.compare %>%
         filter(Year=="2018" & ifelse(Country=="Germany",Age<75,Age<80)),
       aes(Age,Value, group=Type, color=Type,shape=Type, fill=Type))+
  geom_line(aes( linetype=Type,color=Type), size=0.8)+
  geom_jitter(aes(shape=Type, color=Type, size=Type))+
  scale_shape_manual(values=c(1,21, 17, 17, 17,17))+
  scale_color_manual(values=c("black",'blue','brown', 'black',"green","pink"))+
  scale_size_manual(values=c(3,1,1,1,1,1))+
  #scale_color_viridis_d()+
  scale_linetype_manual(values=c("blank","blank", "solid","solid","solid","solid"))+
  facet_grid(Sex~Country)+
  geom_rangeframe() +
  #theme_tufte(base_size = 20)+
  theme_bw(base_size = 20)+
  theme(legend.position="bottom")+ #c(0.9,0.2),
  #strip.text.x = element_blank())+
  ylab("Age-Specific Prevalence")


# heat map to incorporate all years and smoothing and open-end interval at 80+

X11(width=30, height=25)

pal <- wes_palette("Zissou1", 103, type = "continuous")

ggplot(prev.compare %>%
         filter(Age>20& Country!="Germany"),
       aes(Age,Year,fill=Value))+
  geom_tile(aes(fill=Value),size=100)+
  #scale_fill_distiller(palette = "YlGnBu") +
  #scale_fill_viridis_c(option = "B", direction = -1) +
  scale_fill_gradientn(colours = pal)+
  theme_classic()+
  facet_grid(Country~Type)


# only recent year and breaking by sex
ggplot(prev.compare %>%
         filter(Age>17& Country!="DE", Year%in% c(2018)),
       aes(Age,Country,fill=Value))+
  geom_raster(aes(fill=Value), hjust = 0, vjust = 0)+
  #scale_fill_distiller(palette = "YlGnBu") +
  scale_fill_viridis_c(option = "B", direction = -1) +
  theme_classic()+
  facet_grid(.~Type)


library(gridExtra)

# comparing hly by smoothing method

hly.compare<-HLY.1 %>%
  group_by(Country,Year,Sex,Age) %>%
  pivot_longer(starts_with("HLY.1"), names_to="HLY_type",values_to="HLY",
               names_repair = "unique")

# create long dataset to perform graphs for single age ending in 80+

prev.compare_lx<-HLY.1 %>%
  group_by(Country,Year,Sex,Age) %>%
  pivot_longer(starts_with("Lx"), names_to="Type",values_to="Value",
               names_repair = "unique")

# Comparing the person-years lived total and Sullivan
X11()
prev.compare_lx$Country<-as.factor(prev.compare_lx$Country)

prev.compare_lx$Country<- factor(prev.compare_lx$Country,
                                 levels=c("DE", "DK","ES","FR", "IT","PL","GB"),
                                 labels = c("Germany","Denmark","Spain","France", "Italy","Poland","UK"))

prev.compare_lx$Sex<-as.factor(prev.compare_lx$Sex)

prev.compare_lx$Sex<- factor(prev.compare_lx$Sex, labels = c("Women", "Men"))


ggplot()+
  geom_line(data=prev.compare_lx %>% filter(Type%in%c("Lx","Lx.healthy_half")
                                            & Year==2017 & ifelse(Country=="Germany",Age<75,Age<80)),
            aes(Age,Value,group=factor(Type)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Country)+
  #scale_color_viridis_d(drop = TRUE)+
  geom_area(data=prev.compare_lx %>% filter(Type%in%c("Lx","Lx.healthy_half") &
                                              Year==2017 & ifelse(Country=="Germany",Age<75,Age<80)), aes(Age,Value,
                                                                                                          group=factor(Type), fill=factor(Type)),position = "identity") +
  scale_fill_viridis_d(name="Person-Years Lived", labels=c("Unhealthy","Healthy"), alpha=0.4)+
  # scale_fill_manual(values=alpha(c('#882255','#009988'),0.5), name="Person-Years Lived", labels=c("Unhealthy","Healthy"))+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y="Person-Years Lived (Lx)", x="Age")


# Comparing the person-years lived for all methods

ggplot()+
  geom_line(data=prev.compare_lx %>% filter(Type%in%c("Lx","Lx.healthy_half","Lx.healthy.spline_half",
                                                      "Lx.healthy.penalty","Lx.healthy.poly")
                                            & Year==2017 & ifelse(Country=="Germany",Age<75,Age<80)),
            aes(Age,Value,group=factor(Type)), color="black",linetype="dashed", size=1.6)+
  facet_grid(Sex~Country)+
  scale_color_viridis_d(drop = TRUE, guide="none")+
  geom_area(data=prev.compare_lx %>% filter(Type%in%c("Lx","Lx.healthy_half","Lx.healthy.spline_half",
                                                      "Lx.healthy.penalty","Lx.healthy.poly")
                                            &  Year==2017& ifelse(Country=="Germany",Age<75,Age<80)), aes(Age,Value,
                                                                                                          group=factor(Type),color=factor(Type), fill=factor(Type)),position = "identity") +
  scale_fill_manual(values = alpha(rev(suf_palette("hanwell", n = 5,type = "continuous")), 0.4),name = "Indicator",
                    labels=c("Unhealthy","Healthy Penalty","Healthy Polynomial","Healthy Splines","Healthy Unsmoothed"))+
  theme_pander(base_size = 14)+
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12))+
  #strip.text.x = element_blank())+
  labs(y="Person-Years Lived (Lx)", x="Age")


# bidimensional smoothing from Camarda for these countries until age 80
# Total Prevalence- doing this for estimating the lists, since Camarda´s package only accepts matrix-like formats
# and a very restricted data format. Check https://rdrr.io/cran/MortalitySmooth/man/MortalitySmooth-package.html

# exposure - total limited + healthy - or the denominator of the prevalence
hly.smooth_exp<-hly.80 %>%
  filter(Country!="DE"& Year%in%c(2008:2018) &Age>=17) %>%
  mutate(exp=Limited+Healthy,
         Year=as.character(Year)) %>%
  arrange(Country,Year,Sex)%>%
  group_by(Country,Sex)%>%
  group_split()%>%
  map( ~ .x %>%
         pivot_wider(1:4,18,names_from=Year,
                     values_from=exp))

# limited - the numerator of prevalence or the number of persons with at least one limitation
hly.smooth_lim<-hly.80 %>%
  filter(Country!="DE"& Year%in%c(2008:2018) &Age>=17) %>%
  mutate(lim=Limited,
         Year=as.character(Year)) %>%
  arrange(Country,Year,Sex)%>%
  group_by(Country,Sex) %>%
  group_split()%>%
  map( ~ .x %>%
         pivot_wider(1:4,18,names_from=Year,
                     values_from=lim ))

# prevalences

hly.smooth_prev<-hly.80 %>%
  filter(Country!="DE"& Year%in%c(2008:2018) &Age>=17) %>%
  mutate(Year=as.character(Year)) %>%
  arrange(Country,Year,Sex)%>%
  group_by(Country,Sex)%>%
  group_split()%>%
  map( ~ .x %>%
         pivot_wider(1:4,18,names_from=Year,
                     values_from=Prev.1))

# extract values for some countries (bigger ones), only women first
# France
fra.exp_80<-hly.smooth_prev[[5]]
names(fra.exp_80)[4:14] <- c(2008:2018)
fra.exp_80<-fra.exp_80 %>%
  select(-c(1:3))
fra.exp_80<-as.matrix(fra.exp_80)
fra.lim_80<-hly.smooth_lim[[5]]
names(fra.lim_80)[4:14] <- c(2008:2018)
fra.lim_80<-as.matrix(fra.lim_80)

# Poland
pol.exp_80<-hly.smooth_exp[[9]]
names(pol.exp_80)[4:14] <- c(2008:2018)
pol.lim_80<-hly.smooth_lim[[9]]
names(pol.lim_80)[4:14] <- c(2008:2018)

# Spain
spain.exp_80<-hly.smooth_exp[[3]]
names(spain.exp_80)[4:14] <- c(2008:2018)
spain.lim_80<-hly.smooth_lim[[3]]
names(spain.lim_80)[4:14] <- c(2008:2018)
# UK
uk.exp_80<-hly.smooth_exp[[11]]
names(uk.exp_80)[4:14] <- c(2008:2018)
uk.lim_80<-hly.smooth_lim[[11]]
names(uk.lim_80)[4:14] <- c(2008:2018)

uk.exp_80<-uk.exp_80 %>%
  select(-c(1:3))
uk.exp_80<-as.matrix(uk.exp_80)
uk.lim_80<-hly.smooth_lim[[11]]
names(uk.lim_80)[4:14] <- c(2008:2018)
uk.lim_80<-as.matrix(uk.lim_80)


Age<-17:80
Year<-2008:2018

#France
fit2D <- Mort2Dsmooth(x = Age, y = Year, Z = fra.exp_80)
f.plot<-plot(fit2D, palette = "terrain.colors")

X11()
f.plot

#UK
fit2D_uk <- Mort2Dsmooth(x = Age, y = Year, Z = uk.exp_80)
f.plot_uk<-plot(fit2D, palette = "terrain.colors")

X11()
f.plot_uk


## since the aim is to interpolate FEW data-points
## we use a large number of B-splines allows a precise,
## but not parsimonius, description

# average changes
#HT_fra_growth<-HT_fra %>%
#  filter(Age %in% c(0,50,65) & Year>2007) %>%
#  select(1:4,20,21) %>%
#  group_by(Country, Age, Sex, HLY.Indicator) %>%
#  mutate(change=(HLY-lag(HLY))/lag(HLY)*100) %>%
#  drop_na() %>%
#  ungroup() %>%
#  group_by(Country, Age, Sex, HLY.Indicator) %>%
#  mutate(avg.change=mean(change))
