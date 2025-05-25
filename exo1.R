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
  geom_boxplot(aes(x=variable, y=log(value), color=variable)) +
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

# Figure la plus utile !!!
library(GGally)
ggpairs(raisin[,-p])
ggpairs(raisin, ggplot2::aes(colour=as.factor(Class), alpha=0.4))

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
valp
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
# stopifnot(max(abs(F_manual[,1:7] - pca_all$ind$coord)) < 1e-8) # identique

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
                                
head(round(contrib_ind, 3))

# --- Méthode 2 : via FactoMineR + factoextra pour vérification et visualisation
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



# FAIRE cos2 pour variables
# FAIRE représentativité des variables et individus



## 3. CAH Ward (style TP5) ----------------------------------------------

library(cluster)      # silhouette
library(factoextra)   # fviz_dend, fviz_cluster, fviz_silhouette
library(alluvial)     # alluvial
library(mclust)       # adjustedRandIndex

# a) Préparation et CAH
X.sc    <- scale(raisin[,-8]) 
res.ward <- hclust(dist(X.sc), method="ward.D2")

# b) Dendrogramme + éboulis des hauteurs
par(mfrow = c(1,2), mai = c(0.6,0.6,0.6,0.2))
# dendrogramme avec rectangles
plot(res.ward, labels = FALSE, hang = -1)
rect.hclust(res.ward, k = 2, border = pal)

# éboulis des hauteurs de fusion
barplot(rev(res.ward$height)[1:20], 
        main = "Éboulis des hauteurs (Ward)", 
        xlab = "fusions récentes", 
        ylab = "hauteur")
abline(h = rev(res.ward$height)[2], lty = 2) 
mtext("seuil k=2", side=3, line=-2, at=1)

# c) silhouette moyenne pour k = 2:6
sil_vals <- sapply(2:6, function(k){
  mean(silhouette(cutree(res.ward,k), dist(X.sc))[,3])
})
plot(2:6, sil_vals, type="b", pch=19, 
     xlab="Nombre de groupes k", ylab="Silhouette moyenne",
     main="Évolution de la silhouette (Ward)")
k_opt <- which.max(sil_vals) + 1
abline(v = k_opt, col="red", lty=2)
text(k_opt, max(sil_vals), paste0("k*=",k_opt), pos=3, col="red")

# d) Analyse en k = 2 ----------------------------------------------------
cl2 <- cutree(res.ward, 2)

# erreur de classification et ARI
err2 <- mean(cl2 != as.numeric(raisin$Class))
ARI2 <- adjustedRandIndex(cl2, as.numeric(raisin$Class))
cat("k=2 → erreur =",round(err2,3),
    "| ARI =",round(ARI2,3), "\n")

# silhouette détaillée
sil2 <- silhouette(cl2, dist(X.sc))
par(mfrow=c(1,1))
fviz_silhouette(sil2) + 
  ggtitle(paste("Silhouette CAH (Ward, k=2) – avg =", round(mean(sil2[,3]),3)))

# indice de Calinski–Harabasz
inter_inert <- function(x, g){
  ni <- table(g); mu <- colMeans(x)
  sum(sapply(names(ni), function(cl){
    nci <- ni[cl]
    mui <- colMeans(x[g==cl,,drop=FALSE])
    nci * sum((mui - mu)^2)
  })) / nrow(x)
}
B2   <- inter_inert(X.sc, cl2) * nrow(X.sc)
Wtot <- sum(rowSums(X.sc^2))
W2   <- Wtot - B2
CH2  <- (B2/(2-1)) / (W2/(nrow(X.sc)-2))
cat("Calinski–Harabasz (k=2) =", round(CH2,3), "\n")

# e) Flux alluvial vs. classe vraie
library(alluvial)
# on prépare les deux partitions : vraie classe vs CAH k=2
df_alluv <- data.frame(
  TrueClass = factor(raisin$Class),
  Ward2     = factor(cl2, labels = c("1","2"))
)

# on dessine sans blocs, en utilisant pal pour les couleurs
alluvial(
  df_alluv,
  freq   = rep(1, nrow(df_alluv)),
  col    = pal[as.numeric(df_alluv$TrueClass)],
  border = pal[as.numeric(df_alluv$TrueClass)],
  hide   = FALSE,
  alpha  = 0.6)


## 4. Clustering sur k premières CP --------------------------------------------------------------
library(FactoMineR)   # PCA
library(mclust)       # adjustedRandIndex

# 1) ACP FactoMineR (centrage-réduction intégrée)
pca_res <- PCA(raisin[ , -8], scale.unit = TRUE, ncp = 7, graph = FALSE)

# 2) Coordonnées des individus
coords <- pca_res$ind$coord    # matrice n × 7

# 3) Boucle sur k = 1…7 pour CAH Ward + évaluation
results <- t(sapply(1:7, function(k) {
  hc  <- hclust(dist(coords[, 1:k]), method = "ward.D2")
  cl  <- cutree(hc, k = 2)
  err_raw <- mean(cl != as.numeric(raisin$Class))
  err     <- min(err_raw, 1 - err_raw)             # erreur corrigée
  ari     <- adjustedRandIndex(cl, as.numeric(raisin$Class))
  c(err_raw = err_raw, err = err, ARI = ari)
}))
colnames(results) <- c("err_raw", "err_corr", "ARI")
results

# 4) Affichage de l’erreur corrigée et de l’ARI
par(mfrow = c(1,2))
plot(1:7, results[,"err_corr"], type = "b", pch = 19, col = "steelblue",
     xlab = "k (CP)", ylab = "Erreur corrigée",
     main = "Erreur CAH vs. k premières CP")
k_opt <- which.min(results[,"err_corr"])
abline(v = k_opt, col = "red", lty = 2)
text(k_opt, results[k_opt,"err_corr"], paste0("k*=", k_opt), pos = 3, col = "red")

plot(1:7, results[,"ARI"], type = "b", pch = 19, col = "darkgreen",
     xlab = "k (CP)", ylab = "ARI",
     main = "Adjusted Rand Index vs. k")

cat(sprintf("k optimal = %d  (erreur corrigée = %.3f, ARI = %.3f)\n",
            k_opt,
            results[k_opt,"err_corr"],
            results[k_opt,"ARI"]))

