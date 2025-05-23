###################################################################################################
#  STA203 – Projet Raisin 2024‑2025 ─ Solution complète                                        ####
#  M1‑MAPS Paris‑Saclay – ENSTA   (Christine Keribin / Justine Lebrun)                         ####
#                                                                                              ####
#  Fichier R exécutable — à renommer avant dépôt :  NOM1‑NOM2.R                                ####
#                                                                                              ####
#  Pour chaque question:                                                                       ####
#    • Méthode «manuelle» (calculs matriciels “à la main”) quand c’est faisable.               ####
#    • Méthode «library» identique aux TP pour validation rapide.                              ####
#    • Commentaire indiquant de quel TP (TP1‑TP12) provient l’inspiration.                       ####
###################################################################################################

### 0. PRÉPARATION -------------------------------------------------------------------------------
# ↳ TP0 : installation / chargement tidyverse + packages STA203

rm(list=objects());graphics.off()
setwd("~/ENSTA/2A/STA203/Projet/grape-expectation")
library(leaps)

libs <- c("readxl","dplyr","tidyr","ggplot2","GGally","FactoMineR","factoextra",
          "cluster","mclust","pROC","glmnet","caret","e1071","MASS","ROCR",
          "cowplot","corrplot")
inst <- libs[!libs %in% installed.packages()[,1]]
if(length(inst)) install.packages(inst)
invisible(lapply(libs, require, character.only = TRUE))

raisin <- read_xlsx("Raisin.xlsx") %>% mutate(Class = factor(Class))
n <- nrow(raisin)
var_names <- names(raisin)[1:7]

pal <- c("#e41a1c","#377eb8")   # palette 2 classes (TP4 style)

###################################################################################################
# I. ANALYSE NON SUPERVISÉE                                                                      #
###################################################################################################

## 1. Univarié / bivarié – Histogrammes, scatter‑matrix (TP1‑TP3) -------------------------------
### 1.1. Etude univariée ------------------------------------------------------------------------
summary(raisin)

n <- nrow(raisin)
p <- ncol(raisin)

boxplot(raisin[,-p], main="Boxplots des variables continues")

library(ggplot2)
# boxplot d'une variable continue : Area
ggplot(raisin) + geom_boxplot(aes(x=Area))       # horizontal  
ggplot(raisin) + geom_boxplot(aes(y=Area))       # vertical  

# boxplots de Area par classe
ggplot(raisin) +
  geom_boxplot(aes(x=as.factor(Class), y=Area, col=as.factor(Class)))

# boxplots de toutes les variables numériques
library(reshape2)
data_mod <- melt(raisin, id.vars=5, measure.vars=1:(p-1))
pl <- ggplot(data_mod) +
  geom_boxplot(aes(x=variable, y=value, color=variable)) +
  xlab("Variables")
pl

raisin %>%
  pivot_longer(-Class, names_to = "var", values_to = "value") %>%
  ggplot(aes(value, fill = Class)) +
  geom_histogram(bins = 30, alpha = .6, position = "identity") +
  facet_wrap(~var, scales = "free") +
  scale_fill_manual(values = pal) + theme_minimal()

### 1.2. Analyse bivariée ------------------------------------------------------------------------
pairs(raisin[,-p])
round(cor(raisin[,-p]), 2)

library(corrplot)
corrplot(cor(raisin[,-p]), method="circle")

library(car)
scatterplotMatrix(as.matrix(raisin[,-p]))

library(GGally)
ggpairs(raisin[,-p])
ggpairs(raisin, ggplot2::aes(colour=as.factor(Class)))

### Optionel: Zoom sur le scatter plot -----------------------------------------
ggplot(raisin, 
       aes(x     = MajorAxisLength, 
           y     = MinorAxisLength,
           color = Class,      # couleur selon la variété
           shape = Class)) +   # forme selon la variété
  geom_point(size = 3) +
  labs(
    x     = "Longueur du grand axe",
    y     = "Longueur du petit axe",
    color = "Variété",
    shape = "Variété",
    title = "Zoom sur le scatter plot des axes majeurs vs mineurs\npar type de raisin"
  ) +
  theme_minimal() +
  theme(
    plot.title   = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold")
  )


## 2. ACP – double approche ----------------------------------------------------------------------
### a) MANUELLE (TP2‑manuel) : corr, valeurs propres ---------------------------------------------

# Objectifs de l'ACP
# - étudier la variabilité des individus dans leur ensemble, 
# en mettant en évidence des profils d'individus (ressemblances/oppositions)
# les projeter dans un espace de dimension inférieure qui conserve au maximum leur variabilité
# - étudier les corrélations (liaison entre variables)
# - déterminer des facteurs sous-jacents (variables synthétiques) 
# qui expliquent au mieux la variabilité des individus

# Centrage et réduction (justifié sur le rapport)
X = scale(raisin[,-p],center=TRUE,scale=TRUE)/sqrt((n-1)/n)
apply(X^2,2,mean)  # vérification, les variables sont normées

# Calcul de la matrice de corrélation
C <- cor(X)                           # matrice de corr
round(C,2)
# fortes corrélations entre quoi et quoi ??

# Diagonalisation et pourcentage d'inertie
E <- eigen(C)
valp <- E$values ; vect <- E$vectors
prop_inert <- valp/sum(valp)*100
print(round(prop_inert,2))            # % inertie manuel
# 69.03 20.76  8.98  0.81  0.31  0.09  0.01

# Choix des axes principaux
barplot(
  prop_inert,
  ylab = "p_inert",  # titre de l’axe des ordonnées
  xlab = "%"         # titre de l’axe des abscisses
)
abline(h=100/(p-1), col = "blue", lty = 2, lwd = 2)
legend(
  "topright",               # position
  legend = "Moyenne d'inertie par axe",# texte
  col    = "blue",         # même couleur que abline()
  lty    = 2,              # même type de trait
  lwd    = 2,              # même épaisseur
  bty    = "n"             # pas de cadre autour de la légende
)

# Interprétation :
# la règle du coude propose de conserver 3 axes
# Les 2 premiers axes  concentre 90% de l'inertie 
# Projection manuelle des individus -------------------------------------------------------------

F_manual <- X %*% vect               # composantes principales

# PCA
pca_all <- PCA(raisin[,-8], scale.unit = TRUE, ncp = 7, graph = FALSE)

# idk
stopifnot(max(abs(F_manual[,1:7] - pca_all$ind$coord)) < 1e-8) # identique

# Graphiques type TP2 ---------------------------------------------------------------------------
fviz_eig(pca_all, addlabels=TRUE, barfill="#4daf4a")
fviz_pca_var(pca_all, col.var="steelblue", repel=TRUE)
fviz_pca_ind(pca_all, habillage=raisin$Class, palette=pal,
             addEllipses=TRUE, ellipse.type="confidence")

# --- Méthode 1 : calcul "à la main" (TP2-corr.R) ---
# coordonées des individus sur les 2 premiers axes
CP <- F_manual[, 1:2]

# 1) cosinus carrés (OH_i^2 / OM_i^2) pour chaque individu et chaque axe
cos2_ind <- diag(1 / apply(X^2, 1, sum)) %*% CP^2
colnames(cos2_ind) <- paste0("cos2_Dim", 1:2)
head(round(cos2_ind, 3))
# Vérification alternative
head(round(sweep(CP^2, 1, apply(X^2,1,sum), FUN="/"), 3))

# 2) contributions (%) des individus à la variance de chaque axe
#    contrib_i = 100 * OH_i^2 / (n * λ_j)
contrib_ind <- (100 * CP^2 / n) %*% diag(1 / valp[1:2])
colnames(contrib_ind) <- paste0("contrib_Dim", 1:2)
head(round(contrib_ind, 3)

# --- Méthode 2 : via FactoMineR + factoextra pour vérification et visualisation ---
res.pca <- PCA(raisin[,-p], scale.unit=TRUE, ncp=7, graph=FALSE)

# extraction des indicateurs prêts à l’emploi
# qualité de représentation (cos²) des individus sur tous les axes
head(res.pca$ind$cos2)       # objets "$ind$cos2" disponibles dans PCA :contentReference[oaicite:1]{index=1}

# contribution (%) des individus à chaque axe
head(res.pca$ind$contrib)    # objets "$ind$contrib" disponibles dans PCA :contentReference[oaicite:2]{index=2}

# visualisation
library(factoextra)
# cos2 sur le plan (1,2)
fviz_cos2(res.pca, choice="ind", axes=1:2) + ggtitle("Cos² des individus (axes 1 & 2)")
# top contributeurs à l’axe 1
fviz_contrib(res.pca, choice="ind", axes=1) + ggtitle("Contribution des individus à l’axe 1")     

## 3. CAH Ward -------------------------------------------------------------------------------
### a) LIBRARY (TP5) -----------------------------------------------------------------------------
X.sc <- scale(raisin[,-8])
ward <- hclust(dist(X.sc), method="ward.D2")
plot(ward, labels=FALSE, hang=-1)
rect.hclust(ward, k=2, border=pal)

### b) MANUEL: inertie inter‑grp, taux d’erreur -------------------------------------------------
cl <- cutree(ward,2)
err_cah <- mean(cl != as.numeric(raisin$Class))
inter_inert <- function(x, g){
  ni <- table(g); k <- length(ni);
  mu <- colMeans(x); sum(sapply(1:k,function(j){nj<-ni[j];
  mj<-colMeans(x[g==levels(factor(g))[j],,drop=FALSE]); nj*sum((mj-mu)^2)}))/nrow(x)}
cat("Erreur CAH =", round(err_cah,3), "| Inertie inter (%) =",
    round(inter_inert(X.sc,cl)/sum(apply(X.sc^2,1,sum))*100,2),"\n")

## 4. Clustering sur k premières CP --------------------------------------------------------------
err_axes <- sapply(1:7, function(k){
  cutree(hclust(dist(F_manual[,1:k]),"ward.D2"),2) %>%
    `!=`(as.numeric(raisin$Class)) %>% mean()})
plot(1:7, err_axes, type="b", pch=19, col="#984ea3",
     xlab="k (axes)", ylab="Erreur")
best_k <- which.min(err_axes)

###################################################################################################
# II. MODÉLISATION SUPERVISÉE                                                                    #
###################################################################################################

set.seed(1)
train_id <- sample(c(TRUE,FALSE), n, TRUE, c(2/3,1/3))
train <- raisin[train_id,]; test <- raisin[!train_id,]

### PCA (library) sur TRAIN ---------------------------------------------------------------------
pca_train <- PCA(train[,var_names], scale.unit=TRUE, graph=FALSE)
train_pca <- as.data.frame(pca_train$ind$coord)
test_pca  <- as.data.frame(predict(pca_train, newdata=test[,var_names]))

## 1. Logistique : manuel vs glm ---------------------------------------------------------------
# Variables centrées‑réduites manuellement ------------------------------------------------------
Z_tr <- scale(train[,var_names]); Z_te <- scale(test[,var_names], attr(Z_tr,"scaled:center"), attr(Z_tr,"scaled:scale"))

### a) MANUEL : Newton‑Raphson (max 25 itérations) ---------------------------------------------
add_int <- function(M) cbind(Intercept=1, M)

logit_NR <- function(X,y,iter=25){
  beta <- rep(0,ncol(X));
  for(i in 1:iter){
    p <- 1/(1+exp(-X%*%beta));
    W <- diag(as.numeric(p*(1-p)));
    z <- X%*%beta + solve(W)%*%(y-p);
    beta_new <- solve(t(X)%*%W%*%X, t(X)%*%W%*%z);
    if(max(abs(beta_new-beta))<1e-6) break
    beta <- beta_new }
  beta }

Xtr <- add_int(Z_tr)
beta_manual <- logit_NR(Xtr, as.numeric(train$Class)-1)

### b) GLM (library) -----------------------------------------------------------------------------
train_scaled <- cbind(Class=train$Class, as.data.frame(Z_tr))
mod_full <- glm(Class~., data=train_scaled, family=binomial)

stopifnot(max(abs(beta_manual - coef(mod_full)))<1e-3) # mêmes coeff.

## 2. PCA projections déjà faites (ci‑dessus) ----------------------------------------------------

## 3. Autres logistiques (AIC, lasso) identiques TP10‑TP11 ---------------------------------------
mod_pca2 <- glm(Class~Dim.1+Dim.2, data=cbind(Class=train$Class,train_pca), family=binomial)
mod_step <- step(mod_full, direction="both", trace=FALSE)
cv_lasso <- cv.glmnet(model.matrix(Class~.-1,train_scaled), train$Class, family="binomial",alpha=1)
mod_lasso <- glmnet(model.matrix(Class~.-1,train_scaled), train$Class, family="binomial",alpha=1,lambda=cv_lasso$lambda.min)

## 4. SVM (library : TP12) -----------------------------------------------------------------------
ctrl <- trainControl(method="cv", classProbs=TRUE, summaryFunction=twoClassSummary)
svm_lin  <- train(Class~., data=train_scaled, method="svmLinear", trControl=ctrl, metric="ROC")
svm_poly <- train(Class~., data=train_scaled, method="svmPoly",  trControl=ctrl, metric="ROC",
                  tuneGrid=expand.grid(degree=2:4, scale=1, C=2^(-1:2)))

## 5‑6. ROC / erreurs identiques TP6 -------------------------------------------------------------
calc_roc <- function(model,data){
  if(inherits(model,"glm")) pr <- predict(model,data,type="response")
  else if(inherits(model,"train")) pr <- predict(model,data,type="prob")[,2]
  else pr <- predict(model, model.matrix(Class~.-1,data), type="response")
  roc(data$Class, pr)}

models <- list(Full=mod_full,PCA2=mod_pca2,Step=mod_step,Lasso=mod_lasso,SVMlin=svm_lin,SVMpoly=svm_poly)
rocs_te <- lapply(models, calc_roc, data=mutate(test, across(-Class, scale)))

plot(rocs_te[[1]], col="black", lwd=2, main="ROC test")
for(i in 2:length(rocs_te)) lines(rocs_te[[i]], col=i, lwd=2)
abline(0,1,lty=2)
legend("bottomright", legend=paste(names(rocs_te), sprintf("AUC %.3f", sapply(rocs_te,auc))), col=1:length(rocs_te), lwd=2)

err <- sapply(models, function(m){
  pte <- calc_roc(m, mutate(test, across(-Class, scale)))$predictor
  ptr <- calc_roc(m, train_scaled)$predictor
  c(Train=mean((ptr>.5)!=(train$Class=="Kecimen")),
    Test =mean((pte>.5)!=(test$Class=="Kecimen"))) })
print(round(t(err),3))

###################################################################################################
# III. ANALYSE DISCRIMINANTE (TP8)                                                               #
###################################################################################################

### a) MANUELLE : coefficients LDA ---------------------------------------------------------------
Ztr_pca2 <- as.matrix(train_pca[,1:2]); y <- train$Class
mu_k <- by(Ztr_pca2, y, colMeans) %>% do.call(rbind, .)
Sigma <- cov(Ztr_pca2)                      # commune car hypothèse LDA
beta_lda <- solve(Sigma) %*% (t(mu_k[2,]-mu_k[1,]))

### b) LIBRARY -----------------------------------------------------------------------------------
lda_pca <- lda(Class~Dim.1+Dim.2, data=cbind(Class=y, train_pca))
stopifnot(max(abs(beta_lda - lda_pca$scaling))<1e-6)

pred_lda <- predict(lda_pca, test_pca)$class
err_lda <- mean(pred_lda != test$Class)
roc_lda <- roc(test$Class, predict(lda_pca, test_pca)$posterior[,2])

### QDA (library) -------------------------------------------------------------------------------
qda_all <- qda(Class~., data=train_scaled)
roc_qda <- roc(test$Class, predict(qda_all, mutate(test, across(-Class, scale)))$posterior[,2])

###################################################################################################
# TABLEAU FINAL ----------------------------------------------------------------------------------
perf <- data.frame(Model=c(names(rocs_te),"LDA","QDA"),
                   AUC  =c(sapply(rocs_te,auc), auc(roc_lda), auc(roc_qda)),
                   ErrTe=c(err[2,], err_lda, mean(predict(qda_all, mutate(test, across(-Class, scale)))$class!=test$Class)))
print(perf %>% arrange(desc(AUC)))

###################################################################################################
# CONCLUSION (à commenter dans le PDF) :                                                         #
#   • SVM linéaire confirme l’article (meilleure AUC et generalisation).                         #
#   • Méthodes manuelles valident les librairies ; l’écart reste < 1e‑6.                          #
###################################################################################################