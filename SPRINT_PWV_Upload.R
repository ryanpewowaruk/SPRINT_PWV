## Set up - Read in and clean dataa

library(haven)
library(ggplot2)
library(Hmisc)
library(nlme)
library(rms)

setwd("C:/Users/pewow/Documents/SPRINT/Treatment and Stiffness Mechanisms")

# Read in stiffness data
MIND <- read.csv("SPRINT_MIND_Stiffness_Mechanisms.csv")

# Read in demographics and baseline descriptors
DEMO_raw <-  read_sas("baseline.sas7bdat")

DEMO <- data.frame(
  maskid = DEMO_raw$maskid,
  age = DEMO_raw$age,
  female = DEMO_raw$female,
  race = DEMO_raw$race4,
  bmi = DEMO_raw$BMI,
  smoke = DEMO_raw$smoke_3cat,
  cvd = DEMO_raw$sub_cvd,
  ckd = DEMO_raw$sub_ckd,
  fbg = DEMO_raw$GLUR,
  n_BP_med = DEMO_raw$n_agents )

DEMO$race2 <- DEMO$race == "BLACK" # Make 2-category race variable Black and non-Black

MIND$year<- as.numeric( sub('.', '', MIND$Visit) ) # Remove 'T' from visit string

# Sort PWV data into baseline (year 0) and follow up (years 1-3)
baseline <- subset(MIND, MIND$year==0)
followup <- subset(MIND, MIND$year>0)
followup <- data.frame(maskid = followup$maskid,
                      year = followup$year,
                      PWV = followup$PWV,
                      PWV_struct = followup$PWV_struct,
                      PWV_LD = followup$PWV_LD,
                      SP = followup$SP,
                      DP = followup$DP)

# Rename baseline data to be tagged with a "0" at the end
baseline$year <- NULL
baseline$PWV0 <- baseline$PWV
baseline$PWV_struct0 <- baseline$PWV_struct
baseline$PWV_LD0 <- baseline$PWV_LD
baseline$SP0 <- baseline$SP
baseline$DP0 <- baseline$DP
baseline$PWV <- NULL
baseline$PWV_struct <- NULL
baseline$PWV_LD <- NULL
baseline$SP <- NULL
baseline$DP <- NULL

# Merge baseline stiffness and demographic data
baseline <- merge(baseline, DEMO, by="maskid")

# Merge follow up and baseline data
both <- merge(followup, baseline, by = "maskid")

## Calculate stiffness mechanisms

# Log-mean blood pressure
both$LM_BP <- with(both, (SP-DP)/log(SP/DP) )
both$LM_BP0 <- with(both, (SP0-DP0)/log(SP0/DP0) )

# Log-mean BP and 120/80 mmHg
LM_12080 <- (120-80)/log(120/80)

# Total PWV at SBP/DBP
both$PWV_Tot <- with(both, sqrt( (LM_BP/DP) * PWV^2 + 133.32*LM_BP/1050 * log(LM_BP/DP) ) )  
# Structural PWV at 120/80
both$PWV_12080 <- with(both, sqrt( (LM_12080/DP) * PWV^2 + 133.32*LM_12080/1050 * log(LM_12080/DP) ) )  
# Load-Dependent PWV = Total - Structural
both$PWV_12080_LD <- with(both, PWV_Tot - PWV_12080 )

# Same calculations for baseline
both$PWV_Tot0 <- with(both, sqrt( (LM_BP0/DP0) * PWV0^2 + 133.32*LM_BP0/1050 * log(LM_BP0/DP0) ) )  
both$PWV_12080_0 <- with(both, sqrt( (LM_12080/DP0) * PWV0^2 + 133.32*LM_12080/1050 * log(LM_12080/DP0) ) )  
both$PWV_12080_LD0 <- with(both, PWV_Tot0 - PWV_12080_0 )

# Sort baseline to only include participants with follow up data
baseline_w_FU <- subset(baseline, maskid %in% both$maskid)

# Calculate stiffness mechanisms for these participants
baseline_w_FU$LM_BP0 <- with(baseline_w_FU, (SP0-DP0)/log(SP0/DP0) )
baseline_w_FU$PWV_Tot0 <- with(baseline_w_FU, sqrt( (LM_BP0/DP0) * PWV0^2 + 133.32*LM_BP0/1050 * log(LM_BP0/DP0) ) )  
baseline_w_FU$PWV_12080_0 <- with(baseline_w_FU, sqrt( (LM_12080/DP0) * PWV0^2 + 133.32*LM_12080/1050 * log(LM_12080/DP0) ) )  
baseline_w_FU$PWV_12080_LD0 <- with(baseline_w_FU, PWV_Tot0 - PWV_12080_0 )

# Sort follow up data into years 1-3
year1 <- subset(both, year==1)
year2 <- subset(both, year==2)
year3 <- subset(both, year==3)

# Define datadist so RMS functions work
dd <- datadist(both)
options(datadist='dd')


## Make Table 1
# Table 1 variables
Table1 <- data.frame(age = baseline_w_FU$age,
                     female = baseline_w_FU$female,
                     race = baseline_w_FU$race,
                     bmi = baseline_w_FU$bmi,
                     smoke = baseline_w_FU$smoke,
                     SP = baseline_w_FU$SP0,
                     DP = baseline_w_FU$DP0,
                     n_BP_med = baseline_w_FU$n_BP_med,
                     cvd = baseline_w_FU$cvd,
                     ckd = baseline_w_FU$ckd,
                     Tx = baseline_w_FU$randAssign)

# Standard and intensive treatment groups
Table1_Std <- subset(Table1, Tx==0)
Table1_Int <- subset(Table1, Tx==1)

describe(Table1_Std)                  
describe(Table1_Int)                  

## Table 2: Descriptive Longitudinal PWV Data

# Baseline
describe( subset(baseline_w_FU, randAssign==0, select=c(SP0, DP0, PWV0, PWV_Tot0, PWV_12080_0, PWV_12080_LD0) ) )
describe( subset(baseline_w_FU, randAssign==1, select=c(SP0, DP0, PWV0, PWV_Tot0, PWV_12080_0, PWV_12080_LD0) ) )

# Year 1
describe( subset(year1, randAssign==0, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )
describe( subset(year1, randAssign==1, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )

# Year 2
describe( subset(year2, randAssign==0, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )
describe( subset(year2, randAssign==1, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )

# Year 3
describe( subset(year3, randAssign==0, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )
describe( subset(year3, randAssign==1, select=c(SP, DP, PWV, PWV_Tot, PWV_12080, PWV_12080_LD) ) )


## Determine correlation structure

cp <- list(corCAR1,corCompSymm,corGaus)
test_corr  <- vector('list',length(cp))

# Systolic BP
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(SP ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(SP0, 4), data=both,
                        correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))

# Diastolic BP
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(DP ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(DP0, 4), data=both,
                        correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))

# Measured cfPWV
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(PWV ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV0, 4), data=both,
                                correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))


# Total PWV
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(PWV_Tot ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_Tot0, 4), data=both,
                        correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))

# Structural PWV
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(PWV_12080 ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                        correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))

# Load-Dependent PWV
for (j in 1:length(cp)) {
  
  test_corr[[j]] <- Gls(PWV_12080_LD ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                        correlation=cp[[j]](form=~year | maskid))
}
c(AIC(test_corr[[1]]), AIC(test_corr[[2]]), AIC(test_corr[[3]]))
c(BIC(test_corr[[1]]), BIC(test_corr[[2]]), BIC(test_corr[[3]]))

# Compount symmetry had best fit correlation structure for all outcome variables

## Longitudinal Model for SBP
a_SP <- Gls(SP ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(SP0, 4), data=both,
               correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_SP); both$fitted <- fitted(a_SP)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=SP0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

# Plot spline over time
ylm <- coord_cartesian(ylim = c(120, 150))
ggplot(Predict(a_SP, year=seq(1, 3, by=0.1), randAssign), aes(size=2),
       adj.subtitle=FALSE, legend.position='top') + ylm +
  geom_line(size = 1)+
  annotate("line", x=c(0, 1), y=c(139, 139.4), color="black", linetype=2, linewidth=1 ) +
  annotate("line", x=c(0, 1), y=c(139, 127.6), color="red" , linetype=2, linewidth=1 ) +
  labs(title = "A. Systolic BP", x = "Year", y = "BP (mmHg)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw() +
  theme( text = element_text(size = 18) )

# Calculate treatment contrasts for years 1-3
k_SP <- contrast(a_SP, list(year=c(1, 2, 3), randAssign=1),
                    list(year=c(1, 2, 3), randAssign=0))
print(k_SP, digits=3)

## Longitudinal Model for DBP
a_DP <- Gls(DP ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(DP0, 4), data=both,
            correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_DP); both$fitted <- fitted(a_DP)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=DP0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

## Plot spline over time
ggplot(Predict(a_DP, year=seq(1, 3, by=0.1), randAssign),
       adj.subtitle=FALSE, legend.position='top') +
  annotate("line", x=c(0, 1), y=c(77, 76.15), color="black", linetype=2 ) +
  annotate("line", x=c(0, 1), y=c(77, 70.53), color="red" , linetype=2 ) +
  labs(title = "Diastolic BP", x = "Year", y = "DBP (mmHg)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw()

# Calculate treatment contrasts for years 1-3
k_DP <- contrast(a_DP, list(year=c(1, 2, 3), randAssign=1),
                 list(year=c(1, 2, 3), randAssign=0))
print(k_DP, digits=3)

## Longitudinal Model for Measured cfPWV
a_meas <- Gls(PWV ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV0, 4), data=both,
         correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_meas); both$fitted <- fitted(a_meas)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=PWV0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

## Plot spline over time
ylm <- coord_cartesian(ylim = c(9, 12))
ggplot(Predict(a_meas, year=seq(1, 3, by=0.1), randAssign),
       adj.subtitle=FALSE, legend.position='top') + ylm +
  annotate("line", x=c(0, 1), y=c(10.4, 10.015), color="black", linetype=2 ) +
  annotate("line", x=c(0, 1), y=c(10.4, 9.826), color="red" , linetype=2 ) +
labs(title = "Measured cfPWV", x = "Year", y = "PWV (m/s)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw()

# Calculate treatment contrasts for years 1-3
k_meas <- contrast(a_meas, list(year=c(1, 2, 3), randAssign=1),
               list(year=c(1, 2, 3), randAssign=0))
print(k_meas, digits=3)

## Longitudinal Model for Total cfPWV @ SBP/DBP
a_total <- Gls(PWV_Tot ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_Tot0, 4), data=both,
               correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_total); both$fitted <- fitted(a_total)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=PWV_Tot0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

## Plot spline over time
ylm <- coord_cartesian(ylim = c(11, 14))
ggplot(Predict(a_total, year=seq(1, 3, by=0.1), randAssign), aes(size=2),
       adj.subtitle=FALSE) + ylm +
  geom_line(size = 1)+
  annotate("line", x=c(0, 1), y=c(12.34, 12.008), color="black", linetype=2, linewidth=1 ) +
  annotate("line", x=c(0, 1), y=c(12.34, 11.756), color="red" , linetype=2, linewidth=1 ) +
  labs(title = "B. Total PWV", x = "Year", y = "PWV (m/s)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw() + theme( text = element_text(size = 18), legend.position='right' )

# Calculate treatment contrasts for years 1-3
k_total <- contrast(a_total, list(year=c(1, 2, 3), randAssign=1),
                    list(year=c(1, 2, 3), randAssign=0))
print(k_total, digits=3)

## Longitudinal Model for Structural PWV 12080
a_struct <- Gls(PWV_12080 ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4) , data=both,
                correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_struct); both$fitted <- fitted(a_struct)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=PWV_12080_0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

## Plot spline over time
ylm <- coord_cartesian(ylim = c(11, 14))
ggplot(Predict(a_struct, year=seq(1, 3, by=0.1), randAssign), aes(size=2),
       adj.subtitle=FALSE, legend.position='right') + ylm +
  geom_line(size = 1)+
  annotate("line", x=c(0, 1), y=c(11.94, 11.70), color="black", linetype=2, linewidth=1 ) +
  annotate("line", x=c(0, 1), y=c(11.94, 12.01), color="red" , linetype=2, linewidth=1 ) +
  labs(title = "C. Structural PWV", x = "Year", y = "PWV (m/s)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw() + theme( text = element_text(size = 18) )


# Calculate treatment contrasts for years 1-3
k_struct <- contrast(a_struct, list(year=c(1, 2, 3), randAssign=1),
                     list(year=c(1, 2, 3), randAssign=0))
print(k_struct, digits=3)


## Longitudinal Model for Load-Dep PWV @ 120/80
a_LD <- Gls(PWV_12080_LD ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
            correlation=corCompSymm(form=~year | maskid))

# Check residuals visually
both$resid <- r <- resid(a_LD); both$fitted <- fitted(a_LD)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
  facet_grid(~ randAssign) + yl
p2 <- ggplot(both, aes(x=PWV_12080_LD0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=year, y=resid)) + yl + ylim(-20,20) +
  stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
  geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)

## Plot spline over time
ylm <- coord_cartesian(ylim = c(-1, 1))
ggplot(Predict(a_LD, year=seq(1, 3, by=0.1), randAssign), aes(size=2),
       adj.subtitle=FALSE, legend.position='top') + ylm +
  geom_line(size = 1)+
  annotate("line", x=c(0, 1), y=c(0.395, 0.35), color="black", linetype=2, linewidth=1 ) +
  annotate("line", x=c(0, 1), y=c(0.395, -0.23), color="red" , linetype=2, linewidth=1 ) +
  labs(title = "D. Load-Dependent PWV", x = "Year", y = "PWV (m/s)", color="") +
  scale_color_manual(labels = c("Standard", "Intensive"), values = c("black", "red")) +
  theme_bw() +
  theme( text = element_text(size = 18) )

# Calculate treatment contrasts for years 1-3
k_LD <- contrast(a_LD, list(year=c(1, 2, 3), randAssign=1),
                 list(year=c(1, 2, 3), randAssign=0))
print(k_LD, digits=3)

## Examine if interactions for sub-groups are overfitting the data

## Structural PWV

# Main effect of age
m_struct_age <- Gls(PWV_12080 ~ rcs(age,3) + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                   correlation=corCompSymm(form=~year | maskid))
# Treatment interaction with age
i_struct_age <- Gls(PWV_12080 ~ rcs(age,3) * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                   correlation=corCompSymm(form=~year | maskid))

# Is including treatment interaction overfitting? Assess with AIC and BIc
c(AIC(m_struct_age), AIC(i_struct_age))
c(BIC(m_struct_age), BIC(i_struct_age))

# Calculate contrast
contrast(i_struct_age, list(randAssign=1, age=c(66, 74, 79), year=3), 
         list(randAssign=0, age=c(66, 74, 79), year=3))

contrast(i_struct_age, list(randAssign=1, age=66, year=3), 
         list(randAssign=0, age=66, year=3),
         list(randAssign=1, age=74, year=3),
         list(randAssign=0, age=74, year=3), type='joint')

contrast(i_struct_age, list(randAssign=1, age=66, year=3), 
         list(randAssign=0, age=66, year=3),
         list(randAssign=1, age=79, year=3),
         list(randAssign=0, age=79, year=3), type='joint')

# Repeat for sex, prior CVD, CKD, race, and baseline PWV
m_struct_sex <- Gls(PWV_12080 ~ female + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

# Repeat for sex, prior CVD, CKD, race, and baseline PWV
i_struct_sex <- Gls(PWV_12080 ~ female * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_struct_sex), AIC(i_struct_sex))
c(BIC(m_struct_sex), BIC(i_struct_sex))

contrast(i_struct_sex, list(randAssign=1, female=c(0,1), year=3), 
         list(randAssign=0, female=c(0,1), year=3))

contrast(i_struct_sex, list(randAssign=1, female=0, year=3), 
         list(randAssign=0, female=0, year=3),
         list(randAssign=1, female=1, year=3),
         list(randAssign=0, female=1, year=3), type='joint')

Predict(i_struct_sex, year=3, randAssign, female)

m_struct_cvd <- Gls(PWV_12080 ~ cvd + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_struct_cvd <- Gls(PWV_12080 ~ cvd * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_struct_cvd), AIC(i_struct_cvd))
c(BIC(m_struct_cvd), BIC(i_struct_cvd))

contrast(i_struct_cvd, list(randAssign=1, cvd=c(0,1), year=3), 
         list(randAssign=0, cvd=c(0,1), year=3))

contrast(i_struct_cvd, list(randAssign=1, cvd=0, year=3), 
         list(randAssign=0, cvd=0, year=3),
         list(randAssign=1, cvd=1, year=3),
         list(randAssign=0, cvd=1, year=3), type='joint')

m_struct_ckd <- Gls(PWV_12080 ~ ckd + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_struct_ckd <- Gls(PWV_12080 ~ ckd * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_struct_ckd), AIC(i_struct_ckd))
c(BIC(m_struct_ckd), BIC(i_struct_ckd))

contrast(i_struct_ckd, list(randAssign=1, ckd=c(0,1), year=3), 
         list(randAssign=0, ckd=c(0,1), year=3))

contrast(i_struct_ckd, list(randAssign=1, ckd=0, year=3), 
         list(randAssign=0, ckd=0, year=3),
         list(randAssign=1, ckd=1, year=3),
         list(randAssign=0, ckd=1, year=3), type='joint')

m_struct_race <- Gls(PWV_12080 ~ race2 + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_struct_race <- Gls(PWV_12080 ~ race2 * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_struct_race), AIC(i_struct_race))
c(BIC(m_struct_race), BIC(i_struct_race))

contrast(i_struct_race, list(randAssign=1, race2=c(0,1), year=3), 
         list(randAssign=0, race2=c(0,1), year=3))

contrast(i_struct_race, list(randAssign=1, race2=0, year=3), 
         list(randAssign=0, race2=0, year=3),
         list(randAssign=1, race2=1, year=3),
         list(randAssign=0, race2=1, year=3), type='joint')

m_struct_pwv <- Gls(PWV_12080 ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_0, 4), data=both,
                     correlation=corCompSymm(form=~year | maskid))

i_struct_pwv <- Gls(PWV_12080 ~ randAssign * rcs(year, parms=c(1,2,3)) * rcs(PWV_12080_0, 4), data=both,
                     correlation=corCompSymm(form=~year | maskid))

c(AIC(m_struct_pwv), AIC(i_struct_pwv))
c(BIC(m_struct_pwv), BIC(i_struct_pwv))

contrast(i_struct_pwv, list(randAssign=1, PWV_12080_0=c(10.0, 12.0, 14.4), year=3), 
         list(randAssign=0, PWV_12080_0=c(10.0, 12.0, 14.4), year=3))

contrast(i_struct_pwv, list(randAssign=1, PWV_12080_0=10, year=3), 
         list(randAssign=0, PWV_12080_0=10, year=3),
         list(randAssign=1, PWV_12080_0=12, year=3),
         list(randAssign=0, PWV_12080_0=12, year=3), type='joint')

contrast(i_struct_pwv, list(randAssign=1, PWV_12080_0=10, year=3), 
         list(randAssign=0, PWV_12080_0=10, year=3),
         list(randAssign=1, PWV_12080_0=14.4, year=3),
         list(randAssign=0, PWV_12080_0=14.4, year=3), type='joint')

## LD PWV

# Main effect of age
m_ld_age <- Gls(PWV_12080_LD ~ rcs(age,3) + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))
# Treatment interaction with age
i_ld_age <- Gls(PWV_12080_LD ~ rcs(age,3) * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))
# Is including treatment interaction overfitting? Assess with AIC and BIc
c(AIC(m_ld_age), AIC(i_ld_age))
c(BIC(m_ld_age), BIC(i_ld_age))

contrast(i_ld_age, list(randAssign=1, age=c(66, 74, 79), year=3), 
         list(randAssign=0, age=c(66, 74, 79), year=3))

contrast(i_ld_age, list(randAssign=1, age=66, year=3), 
         list(randAssign=0, age=66, year=3),
         list(randAssign=1, age=74, year=3),
         list(randAssign=0, age=74, year=3), type='joint')

contrast(i_ld_age, list(randAssign=1, age=66, year=3), 
         list(randAssign=0, age=66, year=3),
         list(randAssign=1, age=79, year=3),
         list(randAssign=0, age=79, year=3), type='joint')

m_ld_sex <- Gls(PWV_12080_LD ~ female + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_ld_sex <- Gls(PWV_12080_LD ~ female * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_ld_sex), AIC(i_ld_sex))
c(BIC(m_ld_sex), BIC(i_ld_sex))

contrast(i_ld_sex, list(randAssign=1, female=c(0,1), year=3), 
         list(randAssign=0, female=c(0,1), year=3))

contrast(i_ld_sex, list(randAssign=1, female=0, year=3), 
         list(randAssign=0, female=0, year=3),
         list(randAssign=1, female=1, year=3),
         list(randAssign=0, female=1, year=3), type='joint')

m_ld_cvd <- Gls(PWV_12080_LD ~ cvd + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_ld_cvd <- Gls(PWV_12080_LD ~ cvd * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_ld_cvd), AIC(i_ld_cvd))
c(BIC(m_ld_cvd), BIC(i_ld_cvd))

contrast(i_ld_cvd, list(randAssign=1, cvd=c(0,1), year=3), 
         list(randAssign=0, cvd=c(0,1), year=3))

contrast(i_ld_cvd, list(randAssign=1, cvd=0, year=3), 
         list(randAssign=0, cvd=0, year=3),
         list(randAssign=1, cvd=1, year=3),
         list(randAssign=0, cvd=1, year=3), type='joint')

m_ld_ckd <- Gls(PWV_12080_LD ~ ckd + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

i_ld_ckd <- Gls(PWV_12080_LD ~ ckd * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                    correlation=corCompSymm(form=~year | maskid))

c(AIC(m_ld_ckd), AIC(i_ld_ckd))
c(BIC(m_ld_ckd), BIC(i_ld_ckd))

contrast(i_ld_ckd, list(randAssign=1, ckd=c(0,1), year=3), 
         list(randAssign=0, ckd=c(0,1), year=3))

contrast(i_ld_ckd, list(randAssign=1, ckd=0, year=3), 
         list(randAssign=0, ckd=0, year=3),
         list(randAssign=1, ckd=1, year=3),
         list(randAssign=0, ckd=1, year=3), type='joint')

m_ld_race <- Gls(PWV_12080_LD ~ race2 + randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                     correlation=corCompSymm(form=~year | maskid))

i_ld_race <- Gls(PWV_12080_LD ~ race2 * randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                     correlation=corCompSymm(form=~year | maskid))

c(AIC(m_ld_race), AIC(i_ld_race))
c(BIC(m_ld_race), BIC(i_ld_race))

contrast(i_ld_race, list(randAssign=1, race2=c(0,1), year=3), 
         list(randAssign=0, race2=c(0,1), year=3))

contrast(i_ld_race, list(randAssign=1, race2=0, year=3), 
         list(randAssign=0, race2=0, year=3),
         list(randAssign=1, race2=1, year=3),
         list(randAssign=0, race2=1, year=3), type='joint')

m_ld_pwv <- Gls(PWV_12080_LD ~ randAssign * rcs(year, parms=c(1,2,3)) + rcs(PWV_12080_LD0, 4), data=both,
                 correlation=corCompSymm(form=~year | maskid))

i_ld_pwv <- Gls(PWV_12080_LD ~ randAssign * rcs(year, parms=c(1,2,3)) * rcs(PWV_12080_LD0, 4), data=both,
                 correlation=corCompSymm(form=~year | maskid))

c(AIC(m_ld_pwv), AIC(i_ld_pwv))
c(BIC(m_ld_pwv), BIC(i_ld_pwv))

contrast(i_ld_pwv, list(randAssign=1, PWV_12080_LD0=c(-0.19, 0.38, 0.95), year=3), 
         list(randAssign=0, PWV_12080_LD0=c(-0.19, 0.38, 0.95), year=3))

contrast(i_ld_pwv, list(randAssign=1, PWV_12080_LD0=-0.19, year=3), 
         list(randAssign=0, PWV_12080_LD0=-0.19, year=3),
         list(randAssign=1, PWV_12080_LD0=0.38, year=3),
         list(randAssign=0, PWV_12080_LD0=0.38, year=3), type='joint')

contrast(i_ld_pwv, list(randAssign=1, PWV_12080_LD0=-0.19, year=3), 
         list(randAssign=0, PWV_12080_LD0=-0.19, year=3),
         list(randAssign=1, PWV_12080_LD0=0.95, year=3),
         list(randAssign=0, PWV_12080_LD0=0.95, year=3), type='joint')

# All treatment interactions were overfitting based on AIC and BIC metrics.
# This indicates that any heterogeneous treatment effects in our sample are more
# likely noise than signal.