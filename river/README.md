# udn / river score summary


## 1.	Description of the RIVER model.
RIVER (RNA-informed variant effect on regulation) is a hierarchical Bayesian model to infer rare variants of their regulatory effects. Compared with other variant scoring methods, RIVER takes the advantage of utilizing both genomic information and transcriptome information, the formulation of the model is specified in ref1.

## 2.	Model training and variant score generation.
We used GTEx v7 whole genome sequencing and cross-tissue RNA-seq data as training data for the model. The trained model (with learned parameters) is subsequently applied on UDN data to predict effects of rare variants. The model uses rare variants and the genomic annotations at those variants as predictors, uses RNA status (for this case is outlier status based on total gene expression levels) as the target/response variable. 

## 3.	Definition of rare variants.
Rare variants are defined as those of minor allele frequency (MAF) < 0.01 in 1000 genome project phase III all populations combined. For variants in GTEx we additionally require MAF < 0.01 within GTEx cohort itself and for variants in UDN we additionally require MAF < 0.02 within UDN cohort itself. We consider all rare variants 10kb near genes (10kb before transcription start site until 10kb after transcription end site). Overall, there are a median of 2 rare variants per gene and individual pair for GTEx subjects and UDN subjects. For this analysis, we consider protein-coding and lincRNA genes only. 

## 4.	Feature selection.
We used the following genomic annotations: Ensembl VEP, CADD, DANN, conservation score (Gerp, PhyloP, PhastCons), CpG, GC, chromHMM and Encode chromatin-openness track. We selected those features based on their prior evidence of association with regulatory effects (ref1). Features are aggregated over each gene and individual pair, using either max(), min() for quantitative features, or any() for categorical features. TSS distance is inverse-transformed first by taking 1/(1+TSS). Finally, all values are standardized by mean and variance to be used as predictors in the model. 

## 5.	Expression outliers.
Outliers (the response variable) are defined as those with absolute Z-score > 2. Z-scores are calculated based on total gene expression level RPKM from RNA-seq. In addition, for GTEx training data, gene expression levels are corrected by PEER (ref2) to remove technical artifacts and major common-variant eQTL effects are also removed. Z-scores for GTEx are median over all available tissues (ref1). 

Ref1: The impact of rare variation on gene expression across tissues. X Li, Y Kim, EK Tsang, JR Davis. Nature, 2017
Ref2: Stegle, O., Parts, L., Piipari, M., Winn, J. & Durbin, R. Using probabilistic estimation of expression residuals (PEER) to obtain increased power and interpretability of gene expression analyses. Nat. Protoc. 7, 500–507 (2012).
RIVER: https://github.com/ipw012/RIVER
Feature generation: https://github.com/xinli-git/udn



