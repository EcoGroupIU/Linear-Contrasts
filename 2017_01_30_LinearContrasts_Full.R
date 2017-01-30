# Stats and Snacks - "a priori linear contrasts"/"General Linear Hyp. Testing"
# EcoLunch 1/30/17 
# Led by: Briana K. Whitaker

################################################################################
# Helpful Resources :
# Hothorn et al 2015 'multcomp' package documentation
#   https://cran.r-project.org/web/packages/multcomp/multcomp.pdf
# Hothorn et al 2015 additional 'multcomp' examples
#   https://cran.r-project.org/web/packages/multcomp/vignettes/multcomp-examples.pdf
# Bretz/Hothorn/Westfall 2011 Book 'Multiple Comparisons Using R'
#   http://www.ievbras.ru/ecostat/Kiril/R/Biblio/R_eng/Bretz%20Multiple%20Comparisons.pdf

################################################################################

# clear your working environment
rm(list=ls())

# Load Packages into R
library(lme4) 
library(car)
library(ggplot2)
library(multcomp)


# Set Working Directory
#setwd("/Users/brianakwhitaker/Box Sync/Documents/3-Research Projects/2015-Asteraceae-Feedbacks")
setwd("/Users/brianakwhitaker/Box Sync/Documents/5-Teaching&Mentoring")

# Load dataset
data <- read.csv("./2017_01_30_DataForEcoLunch.csv", header=T)    
dim(data)
names(data)

# make the block vector a factor
data$block <- as.factor(data$block)

# Reorder the source category so that the control treatment is last
data$source <- factor(data$source, levels = c("AN", "CA", "EP", "VM", "cont"))
levels(data$source)

# make an interaction term
data$interaction_var <- interaction(data$species, data$source, sep="-")
levels(data$interaction_var)

# Then reorder these as well
data$interaction_var <- factor(data$interaction_var, 
                          levels = c("AN-AN", "AN-CA", "AN-EP", "AN-VM", "AN-cont",        
                                     "CA-AN", "CA-CA", "CA-EP", "CA-VM", "CA-cont", 
                                     "EP-AN", "EP-CA", "EP-EP", "EP-VM", "EP-cont",  
                                     "VM-AN", "VM-CA", "VM-EP", "VM-VM", "VM-cont"))
levels(data$interaction_var)


################################################################################

# Basic Model
lm1 <- lm(sqrt(totalbm) ~ species * source + block, data=data)  
Anova(lm1)      # Note this is Type 2 ANOVA
summary(lm1)
coef(lm1)

# Define a Random Intercepts Model with block as a random intercept
lme1 <- lmer(sqrt(totalbm) ~ species*source + (1|block), data=data, REML=FALSE) 
Anova(lme1)
summary(lme1)
fixef(lme1)

# Define a Random Intercepts Model, but hard-code the interaction term
# really great explanation for why this is a valid approach
# http://stats.stackexchange.com/questions/5250/multiple-comparisons-on-a-mixed-effects-model
lme.interaction <- lmer(sqrt(totalbm) ~  interaction_var + (1|block), data=data, REML=FALSE) 
Anova(lme.interaction)
fixef(lme.interaction)

anova(lme1,lme.interaction) #same DF, AIC, BIC, logLik

# Rationale for running the model this way (lme.interaction)  :
# the lme.interaction model alone does not test for the interaction of factor 
# variable terms, however what it does do is provide more flexibility for making
# a comparison of the joint effects of specific levels within the treatment 
# groups with the 'multcomp' package.



################################################################################

# Introduction to Specific Contrasts - Special Type = Tukey's Post-Hoc
tukey.source <- glht(lme1, linfct=mcp(source="Tukey"))
# functions
?glht   # linfct    # mcp
# function inputs --> 
#    'lme1' which is the model, and 'source' which is the model term of interest
# What types of models accepted? aov() , lm(), lmer(), glm()
summary(tukey.source, test = adjusted("single-step"))


# test adjusted p-values, many options, default is 'single-step', 
# From the authors themselves: "For the homoscedastic normal linear models, the 
#    [default] functions in the package account for the correlations between 
#    test statistics by using the exact multivariate t-distribution. The 
#    resulting procedures are therefore more powerful than the Bonferroni and 
#    Holm methods..."

# Other Options: none, bonferroni, Shaffer, Westfall, single-step, 
        # holm, hochberg, free, hommel


# other pre-exisiting linfct contrast codings are : 
#   Dunnett, Tukey, Sequen, AVE, Changepoint, Williams, Marcus, McDermott, 
#   UmbrellaWilliams, GrandMean


################################################################################

# Define matrix of YOUR OWN specific contrasts --->
    # every row should sum to ZERO 
    # (unless option 'rhs' specified as something other than default of ZERO)

# first, remember the order of your treatments
table(data$source)


# Make a comparison between All Trts. versus the Control group
Qcontrol <- rbind("Live-Sterile" = c(1/4, 1/4, 1/4, 1/4,  -1 ) )   
    # NOTE: use of rbind() function
Qcontrol

# Run the general linear hypotheses model using multcomp package
Qcontrol.comp <- glht(lme1, linfct=mcp(source=Qcontrol))
summary(Qcontrol.comp, test = adjusted("single-step"))
# gives an error, why? --> the glht function uses default encodings, so when 
# there are interactions or covariates present, our Q-matrix is not fully
# explicit about how to deal with those
# the current default method is to generate comparisons for main effects only, 
# ignoring interaction terms and co-variates.
# NOTE: there are other p-value adjustment options available, in these examples, I make no adjustments


# Retrieve confidence intervals for the estimates of the specific comparison
Qcontrol.CIs <- confint(Qcontrol.comp, level=0.95)   #quantile = 1.96
# ALTERNATIVE : two-sided 95% simultaneous confidence intervals for the Bonferoni Test
#  my.calpha <- qt(1-0.05/2/5, 124)  
    # where N=129 minus 5 levels of 'source' = 124, and 5 is for the levels of source
#Qcontrol.CIs <- confint(Qcontrol.comp, calpha=my.calpha) #calpha=2.61
Qcontrol.CIs


# Make a plot of the specific contrasts you estimated
contplot <- qplot(lhs, estimate, data = Qcontrol.CIs, 
      main="95% CI Live-Sterile Contrast", xlab="", 
      ylab ="Estimate", geom = "pointrange", ymin = lwr, ymax = upr) + 
      coord_flip() + geom_hline(yintercept = 0) 
contplot


################################################################################
# more interesting examples

# a useful way to think about your data, get a handle on its structure
table(data$species, data$source)
table(data$interaction_var)


# Define matrix of specific contrasts
QallMat <- rbind(           "AN" = c(1,     -1/3,   -1/3, -1/3, 0,         # For AN
                                    -1/3,    1/3,    0,    0,   0,
                                    -1/3,    0,      1/3,  0,   0,
                                    -1/3,    0,      0,    1/3, 0), #<--NOTE this extra comma, it's important

                            "CA" = c(1/3,   -1/3,    0,    0,   0, 
                                    -1/3,    1,     -1/3, -1/3, 0,          # For CA
                                     0,     -1/3,    1/3,  0,   0,
                                     0,     -1/3,    0,    1/3, 0), #<--also extra comma
                  
                            "EP" = c(1/3,    0,     -1/3,  0,   0,
                                     0,      1/3,   -1/3,  0,   0,
                                    -1/3,   -1/3,    1,   -1/3, 0,          # For EP
                                     0,      0,     -1/3,  1/3, 0), #<--also extra comma
                  
                            "VM" = c(1/3,    0,     0,    -1/3, 0,
                                     0,      1/3,   0,    -1/3, 0,
                                     0,      0,     1/3,  -1/3, 0,          # For VM
                                    -1/3,   -1/3,  -1/3,   1,   0)    )      #NO EXTRA COMMA, but close the rbind parentheses!

# so what does this look like?
View(QallMat)

# Run the general linear hypotheses model using multcomp package
Qall.comp <- glht(lme.interaction, linfct=mcp(interaction_var=QallMat))
summary(Qall.comp, test = adjusted("single-step"))


# Retrieve confidence intervals for the estimates of the specific comparison
Qall.CIs <- confint(Qall.comp, level=0.95)   #quantile = 2.4437
#my.calpha2 <- qt(1-0.05/2/20, 109)  
  # where N=129 minus 20 levels = 109, and 20 is for the levels of interaction_var#
#Qall.CIs <- confint(Qall.comp, calpha=my.calpha2) #calpha=3.10
Qall.CIs


# Make a plot of the specific contrasts you estimated
allplot <- qplot(lhs, estimate, data = Qall.CIs, 
      main="95% CI Home-Away Contrast", xlab="", 
      ylab ="Estimate", geom = "pointrange", ymin = lwr, ymax = upr) + 
      coord_flip() + geom_hline(yintercept = 0) 
allplot


################################################################################


# caveats : spend some time looking at p-values and confidence intervals for YOUR datasets

# end