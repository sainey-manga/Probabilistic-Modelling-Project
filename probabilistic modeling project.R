library('dplyr')
library('gRbase')
library('igraph')
library('gRim')
library('RBGL')
library('pcalg')
library('Rgraphviz')
library('glasso')
library('SIN')
library('mvtnorm') # multivariate Gaussian distribution library('matrixcalc') # operations with matrices
library('mgm')
happiness<- read.csv("/Users/saineymanga/Desktop/2019.csv")
any(is.na(happiness))
attach(happiness)
head(happiness)
clean_happiness<-select(happiness, -c(Overall.rank, Country.or.region))#####droping off the country of origin and the overall rank.
head(clean_happiness)
## CORRELATION MATRIX
library(corrplot)
corr_matrix<-cor(clean_happiness)
corrplot(corr_matrix, type="upper")
corr_matrix
####Data visualization)
library('purrr')
library('tidyr')
library('ggplot2')

clean_happiness %>%
  keep(is.numeric) %>% 
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram()

####Covariance and concentration/partial covariance
S.happiness <- cov.wt(clean_happiness,method="ML")$cov 
K.happiness <- solve(S.happiness) # inverse 
round(K.happiness*100)
PC.happiness <- cov2pcor(S.happiness)
round(PC.happiness*100)
    
####ploting the graphical representation of the model
gen.happiness <- cmod(~Score*GDP.per.capita*Social.support*Perceptions.of.corruption+ 
                      Score*GDP.per.capita*Healthy.life.expectancy+
                      Score*Social.support*Freedom.to.make.life.choices+
                     Generosity*GDP.per.capita*Perceptions.of.corruption+
                     Generosity*Freedom.to.make.life.choices*Perceptions.of.corruption,data=clean_happiness)
gen.happiness
plot(as(gen.happiness,"igraph"))

gen.happiness.graph<- ug(~Score*GDP.per.capita*Social.support*Perceptions.of.corruption+ 
                           Score*GDP.per.capita*Healthy.life.expectancy+
                           Score*Social.support*Freedom.to.make.life.choices+
                           Generosity*GDP.per.capita*Perceptions.of.corruption+
                           Generosity*Freedom.to.make.life.choices*Perceptions.of.corruption)
plot(as(gen.happines.graph,'igraph'))
##In order to make faster the algorithm it is advisable to pass the cliques
cl.happiness <- getCliques(gen.happiness.graph) 
gen.happiness.cl <- cmod(cl.happiness,data=clean_happiness) 
gen.happiness.cl
####Hypothese testing
ciTest_mvn(list(cov=S.happiness,n.obs=nrow(clean_happiness)),#### the independence of happiness given with GDp given other variable
           set=~Score+GDP.per.capita+Social.support+Healthy.life.expectancy+Freedom.to.make.life.choices+Generosity+
             Perceptions.of.corruption)
ciTest_mvn(list(cov=S.happiness,n.obs=nrow(clean_happiness)),#### the independence of happiness given with GDp given other variable
           set=~Score+GDP.per.capita+Social.support+Healthy.life.expectancy+Freedom.to.make.life.choices+Generosity+
             Perceptions.of.corruption, statistic = 'F')

ciTest_mvn(list(cov=S.happiness,n.obs=nrow(clean_happiness)),
           set=~Score+Healthy.life.expectancy+Social.support+GDP.per.capita+Freedom.to.make.life.choices+
             Perceptions.of.corruption+Generosity)
ciTest_mvn(list(cov=S.happiness,n.obs=nrow(clean_happiness)),
           set=~Score+Generosity+Healthy.life.expectancy+Social.support+GDP.per.capita+Freedom.to.make.life.choices+
             Perceptions.of.corruption)
ciTest_mvn(list(cov=S.happiness,n.obs=nrow(clean_happiness)),#### lets remove genorosity
           set=~Score+ Perceptions.of.corruption+Generosity+Healthy.life.expectancy+Social.support+GDP.per.capita+
             Freedom.to.make.life.choices)

###Asymptotic normality of Fisher’s z transform of the partial correlation
cS<-cov2cor(S.happiness)
gaussCItest(1,2,3:6,list(C=cS,n=nrow(clean_happiness)))
gaussCItest(7,2,c(1,3,4,5,6),list(C=cS,n=nrow(clean_happiness)))

#####CONCENTRATION AND REGRESSION
-K.happiness[1,-1]/K.happiness[1,1]
1/K.happiness[1,1] ## residual variance of score
##### MODEL SELECTION#####
sat.happiness <- cmod(~.^., data=clean_happiness)  
aic.happiness<- stepwise(sat.happiness)##base on AIC, note
plot(as(aic.happiness,"graphNEL"),"fdp")
plot(as(aic.happiness, 'igraph'))
bic.happiness<-stepwise(sat.happiness,k=log(nrow(clean_happiness))) ###base on BIC
bic.happiness
plot(as(bic.happiness,"graphNEL"),"fdp")
plot(as(bic.happiness, 'igraph'))
##stepwise method##
test.happiness <- stepwise(sat.happiness, details=1,"test")###Stepwise using the default significance
plot(test.happiness,"neato")
plot(as(test.happiness, "igraph"))###using "igraph" to improve visibility
ind.happiness<-cmod(~.^1,data=clean_happiness)#####another stepwise method using gRim package
set.seed(123)
forw.happiness<-stepwise(ind.happiness,search="headlong",
                        direction="forward",k=log(nrow(clean_happiness)),details=0) 
forw.happiness
plot(forw.happiness,"neato")
plot(as(forw.happiness, 'igraph'))
##convexx optimization method##
res.lasso<-glasso(cS,rho=0.1)
AM <- res.lasso$wi != 0
diag(AM) <- F
g.lasso <- as(AM, "graphNEL")
nodes(g.lasso)<-names(clean_happiness)
glasso.happiness <- cmod(edgeList(g.lasso),data=clean_happiness)
plot(as(glasso.happiness,"igraph"))
##thesholding method
threshold <- .2
Z <- abs(PC.happiness)
Z[Z<threshold] <- 0
diag(Z)<-0
Z[Z>0] <- 1
g.thresh<-as(Z, "graphNEL")
thresh.happiness <- cmod(edgeList(g.thresh), data=clean_happiness) 
thresh.happiness
plot(thresh.happiness, "neato")
plot(as(thresh.happiness, "igraph"))### using 'igraph' to improve the label visibility
##Simultaneous p-values
psin.happiness<-sinUG(S.happiness,n=nrow(clean_happiness)) 
plotUGpvalues(psin.happiness)
gsin.happiness <- as(getgraph(psin.happiness, 0.1), "graphNEL")##thresholding the simultaneous p-values
plot(gsin.happiness, "neato")
plot(as(gsin.happiness, "igraph"))### improving the label visibility
####COMMON EDGES###
commonedges.happiness<-intersection(as(aic.happiness,"graphNEL"), as(bic.happiness,"graphNEL"))
othermodels<-list(test.happiness,forw.happiness, thresh.happiness,gsin.happiness,glasso.happiness)
othermodels<-lapply(othermodels, as, "graphNEL")
for(ii in 1:length(othermodels)){
     commonedges.happiness<-intersection(commonedges.happiness,othermodels[[ii]]) 
}
plot(commonedges.happiness,"fdp")
plot(as(commonedges.happiness, "igraph"))###improving visibility

###CLASSIFICATION ERROR ON THE DATASET###
set.seed(1)
fit_mgm<- mgm(data = clean_happiness, type = c("c",rep("g",13)), levels = c(3,rep(1,13)), k=2, lambdaSel = "CV",
              lambdaFolds= 10, ruleReg="AND", overparameterize = T)
fit_mgm<- mgm(data = clean_happiness, type = rep('g', 7), levels =rep(1, 7), k= 2,  lambdaSel = "CV",
              lambdaFolds= 10, ruleReg="AND", overparameterize = TRUE)
##prediction##
pred_mgm <- predict(fit_mgm, data= clean_happiness, errorCon = c("RMSE","R2"))
pred_mgm$errors
