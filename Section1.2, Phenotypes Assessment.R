library(BiocInstaller)
biocLite()

biocLite("COPDSexualDimorphism.data")   #Install and attach the COPDSexualDimorphism.data package using biocLite. 
library(COPDSexualDimorphism.data)
# add the object expr.meta to your workspace.
data(lgrc.expr.meta)           
names(expr.meta)
head(expr.meta)    
# What is the number of female participants in this study?
table(expr.meta$gender)      # the answer is 110
# What is the median of the distribution of pack years smoked in this cohort (women and men)?
summary(expr.meta$pkyrs)
median(expr.meta$pkyrs)      # the answer is 40
# plots the data based on $pkyrs and $gender
boxplot(pkyrs~gender, data=expr.meta)
hist(expr.meta$pkyrs)
qqnorm(expr.meta$pkyrs)
boxplot(pkyrs~gender, data=expr.meta)

expr.meta$pyp1 = expr.meta$pkyrs+1   #define a positive-valued variable for transformation analysis
library(MASS)
lm1 = lm(pyp1~gender, data=expr.meta)   #tests for a difference in mean pack years (plus 1) between genders
boxcox(lm1)     #plot of the likelihood function for a transformation model
# For what value of lambda does the likelihood reach its highest value for the model lm1?
a<-boxcox(lm1)
a
a$x[which.max(a$y)]       # 0.4646465

#see the effects of the transformation on symmetry and presence of outliers.
lambda= 0.5
boxplot(I(pyp1^lambda)~gender, data=expr.meta)


