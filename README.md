# BUScorrect R package
The BUScorrect R package implements the BUS model to adjust genomic data for batch effects when there are unknown sample subtypes.

## Introduction
High-throughput experimental data are accumulating exponentially in public databases. However, mining valid scientific discoveries from these abundant resources is hampered by technical artifacts and inherent biological heterogeneity. The former are usually termed batch effects, and the latter is often modelled by subtypes. 

Researchers have long been aware that samples generated on different days are not directly comparable. Samples processed at the same time are usually referred to as coming from the same batch. Even when the same biological conditions are measured, data from different batches can present very different patterns. The variation among different batches may be due to changes in laboratory conditions, preparation time, reagent lots, and experimenters [1]. The effects caused by these systematic factors are called batch effects.

Various batch effects correction methods have been proposed when the subtype information for each sample is known [2,3]. Here we adopt a broad definition for subtype. Subtype is defined as a set of samples that share the same underlying genomic profile, in other words biological variability, when measured with no technical artifacts. For instance, groupings such as case and control can be viewed as two subtypes. However, subtype information is usually unknown, and it is often the main interest of the study to learn the subtype for each collected sample, especially in personalized medicine.

Here, the R package BUScorrect fits a Bayesian hierarchical model, the Batch-effects-correction-with-Unknown-Subtypes model (BUS), to correct batch effects in the presence of unknown subtypes [4]. BUS is capable of (a) correcting batch effects explicitly, (b) grouping samples that share similar characteristics into subtypes, (c) identifying features that distinguish subtypes, and (d) enjoying a linear-order computation complexity. After correcting the batch effects with BUS, the corrected value can be used for other analysis as if all samples are measured in a single batch. BUS can integrate batches measured from different platforms and allow subtypes to be measured in some but not all of the batches as long as the experimental design fulfils the conditions listed in [4].

## Installation 
Although this R package is still under review on the Bioconductor platform, you can install its latest version following the R command below.

```
devtools::install_github("XiangyuLuo/BUScorrect")
```

## User's Guide
Please refer to the "BUScorrect\_user\_guide.pdf" vignetee for detailed function instructions.

## References
1. Leek, Jeffrey T., et al. "Tackling the widespread and critical impact of batch effects in high-throughput data." *Nature Reviews Genetics* 11.10 (2010): 733.
2. Johnson, W. Evan, Cheng Li, and Ariel Rabinovic. "Adjusting batch effects in microarray expression data using empirical Bayes methods." *Biostatistics* 8.1 (2007): 118-127.
3. Leek, Jeffrey T., and John D. Storey. "Capturing heterogeneity in gene expression studies by surrogate variable analysis." *PLoS genetics* 3.9 (2007): e161.
4. Luo, Xiangyu, and Yingying Wei. "Batch effects correction with unknown subtypes". *Journal of the American Statistical Association*. Accepted. 
