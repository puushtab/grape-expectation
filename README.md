# grape-expectation

Partie I – Analyse non supervisée
1. Analyses uni- et bi-variée
Objectif : inspecter chaque variable isolément puis explorer les relations deux à deux.

Univariée :

Pour chaque variable quantitative, calculer moyenne, médiane, écart-type (summary(), sd()) et tracer histogrammes et boxplots (cf. TP1 §2–3 TP1-corr).

Pour la variable qualitative (variété), produire une table de comptage et un barplot.

Bivariée :

Nuages de points pour chaque paire de variables avec couleur selon la variété (pairs(), ggpairs()) et matrice de corrélation (corrplot::corrplot()), comme dans TP2 §1 TP2-corr.

Discussion des corrélations fortes (ex. area–perimeter).

2. Méthodologie d’étude de l’ACP
Objectif : suivre rigoureusement les étapes de l’ACP normée.

Centrage et réduction (scale()) TP2-corr.

Calcul de la matrice de corrélation et diagonalisation (eigen()), extraction des valeurs propres et % d’inertie (TP2 §2) TP2-corr.

Choix du nombre d’axes (critère du coude, % cumulé > 80 %).

Calcul des composantes principales :

Coordonnées des individus (F = X %*% vectp) et des variables (cor(X, F)) TP2-corr.

Critères de qualité (cos², contributions) ENSTA-tp2-ACP-suite.

3. Classification hiérarchique ascendante (CAH)
Objectif : segmenter sans supervision, comparer avec les variétés connues.

Lancer CAH sur les 7 variables TP5-corr.

Dendrogramme (plot(hclust(dist(X)))), choix du nombre de clusters (coupure à hauteur, éboulis).

Évaluer la qualité : silhouette (silhouette()) et erreur de classification en comparant cutree() versus la variété ENSTA-tp5-clustering-su….

Tester l’influence de la normalisation (CAH sur les données centrées-réduites vs brutes).

4. CAH sur composantes principales
Objectif : voir l’effet de la réduction de dimension sur la classification.

Pour k = 1,…,7, appliquer CAH aux k premières composantes (programmation R basique).

Calculer l’erreur de classification pour chaque k et identifier le k minimisant l’erreur.

Discussion du biais potentiel (surapprentissage si on choisit k après avoir vu l’erreur).

Partie II – Choix d’un modèle supervisé
1. Régression logistique
Objectif : modéliser la probabilité de variété Besni vs Kecimen.

Définir le modèle logit et rappeler l’impact du centrage-réduction (cf. Cours §8.2).

Ajuster le modèle complet sur les composantes principales extraites (ou directement sur les variables, selon choix) avec glm(..., family=binomial) TP6_cor.

Discussion théorique/pratique sur standardisation des variables (importance pour l’interprétation et convergence).

2. Découpage apprentissage/test
Objectif : évaluer la performance hors échantillon.

Créer l’échantillon d’apprentissage avec set.seed(1); train = sample(...) projetSTA203-2024-2025.

Projeter les données de test sur l’ACP calculée sur l’apprentissage (méthode PCA + predict.PCA).

3. Estimation de modèles
Objectif : comparer plusieurs approches de sélection de variables et de régularisation.

Modèle complet (toutes les composantes ou variables).

Modèle 2 composantes seules.

Sélection pas à pas AIC (stepAIC()) TP10-Choix-C.

Régression Lasso (glmnet(..., alpha=1)) en sélectionnant λ par validation croisée (cv.glmnet()) ENSTA-tp11-Regul.

Choix des hyperparamètres (critère AIC pour stepwise, cv.glmnet pour lasso).

4. Modèles SVM
Objectif : tester des classifieurs alternatifs.

SVM linéaire (svm(..., kernel="linear")) puis SVM polynomial (kernel="polynomial") avec e1071::svm.

Validation croisée pour le coût et le degré (grid search).

5. Courbes ROC et AUC
Objectif : comparer graphiquement et numériquement.

Sur apprentissage et test, tracer les courbes ROC (pROC::roc ou ROCR::performance) TP9.

Superposer la règle aléatoire (ligne diagonale).

Calculer l’AUC pour chaque modèle, l’afficher dans la légende.