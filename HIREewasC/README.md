## A Brief Summary of HIREewas
HIRE is short for HIgh REsolution, as it can detect cell-type-specific risk-CpG sites in epigenome wide association studies (EWAS), while existing EWAS approaches cannot accomplish the "cell-type-specific" function. HIREewas is the accompanied R package to implement the HIRE model. (The corresponding paper link: https://www.biorxiv.org/content/early/2018/09/12/415109 on the bioRxiv)

## Introduction
In epigenome-wide association studies (EWAS), as samples are measured at the bulk level rather than at the single-cell level, 
the obtained methylome for each sample shows the signals aggregated from distinct cell types 
[1-3]. The cellular heterogeneity leads to two main challenges for analyzing EWAS data.

On the one hand, the cell type compositions differ between samples and can be associated with phenotypes 
[2,3]. Both binary phenotypes, such as the diseased or normal status [2], 
and continuous phenotypes, for example, age [3], have been found to affect the cell type compositions. 
As a result, ignoring the cellular heterogeneity in EWAS can lead to a large number of spurious associations 
[3-6]. On the other hand, 
the phenotype may change the methylation level of a CpG site in some but not all of the cell types. 
Identifying the exact cell types that carry the risk-CpG sites can deepen our understandings of disease mechanisms. 
Nevertheless, such identification is challenging because we can only observe the aggregated-level signals.

However, to the best of our knowledge, no existing statistical method for EWAS can detect cell-type-specific associations despite 
the active research on accounting for cell-type heterogeneity. The existing approaches can be categorized into two schools 
[7]: "reference-based" and "reference-free" methods. As the method names indicate, 
the reference-based methods [1,8] require the reference methylation profiles for 
each cell type to be known a priori, while the reference-free methods do not 
depend on any known methylation reference by employing matrix decomposition techniques 
[9] or extracting surrogate variables including principle components as a special case 
[4,6,10,11]. 

Although all of the existing methods aim to address the cellular heterogeneity problem in EWAS and claim whether 
a CpG site is associated with phenotypes at the **aggregated level**, none of them can identify the risk-CpG sites for 
each **individual cell type**, thus missing the opportunity to obtain finer-grained results in EWAS.

We propose a hierarchical model HIRE [12] to identify the association in EWAS at an unprecedented 
HIgh REsolution: detecting whether a CpG site has any associations with the phenotypes in each cell type. 
HIRE not only substantially improves the power of association detection at the aggregated level as compared to the 
existing methods but also enables the detection of risk-CpG sites for individual cell types. HIRE is applicable to EWAS 
with binary phenotypes, continuous phenotypes, or both. 

The user's guide provides step-by-step instructions for the **HIREewas** R package. 
We believe that, by helping biology researchers understand in which cell types the CpG sites are affected by a disease using 
**HIREewas**, HIRE can ultimately facilitate the development of epigenetic therapies by targeting the specifically affected cell types.

## Installation
```
devtools:install_github("XiangyuLuo/HIREewas")
```

## Instructions
Please refer to the HIREewas vignette or the function manual for details.

## References
[1] Eugene Andres Houseman, William P Accomando, Devin C Koestler, Brock C Christensen,
Carmen J Marsit, Heather H Nelson, John K Wiencke, and Karl T Kelsey. DNA methylation
arrays as surrogate measures of cell mixture distribution. BMC Bioinformatics, 13(1):86,
2012.

[2] Yun Liu, Martin J Aryee, Leonid Padyukov, M Daniele Fallin, Espen Hesselberg, Arni
Runarsson, Lovisa Reinius, Nathalie Acevedo, Margaret Taub, Marcus Ronninger, et al.
Epigenome-wide association data implicate DNA methylation as an intermediary of genetic
risk in rheumatoid arthritis. Nature Biotechnology, 31(2):142-147, 2013.

[3] Andrew E Jaffe and Rafael A Irizarry. Accounting for cellular heterogeneity is critical in
epigenome-wide association studies. Genome Biology, 15(2):R31, 2014.

[4] James Zou, Christoph Lippert, David Heckerman, Martin Aryee, and Jennifer Listgarten.
Epigenome-wide association studies without the need for cell-type composition. Nature
Methods, 11(3):309-311, 2014.

[5] Kevin McGregor, Sasha Bernatsky, Ines Colmegna, Marie Hudson, Tomi Pastinen, Aurelie
Labbe, and Celia MT Greenwood. An evaluation of methods correcting for cell-type
heterogeneity in DNA methylation studies. Genome Biology, 17(1):84, 2016.

[6] Elior Rahmani, Noah Zaitlen, Yael Baran, Celeste Eng, Donglei Hu, Joshua Galanter, Sam
Oh, Esteban G Burchard, Eleazar Eskin, James Zou, et al. Sparse PCA corrects for cell
type heterogeneity in epigenome-wide association studies. Nature Methods, 13(5):443-445,
2016.

[7] Andrew E Teschendorff and Caroline L Relton. Statistical and integrative system-level
analysis of DNA methylation data. Nature Reviews Genetics, 2017.

[8] William P Accomando, John K Wiencke, E Andres Houseman, Heather H Nelson, and
Karl T Kelsey. Quantitative reconstruction of leukocyte subsets using DNA methylation.
Genome Biology, 15(3):R50, 2014.

[9] Eugene Andres Houseman, Molly L Kile, David C Christiani, Tan A Ince, Karl T Kelsey,
and Carmen J Marsit. Reference-free deconvolution of DNA methylation data and
mediation by cell composition effects. BMC Bioinformatics, 17(1):259, 2016.

[10] Jeffrey T Leek and John D Storey. Capturing heterogeneity in gene expression studies by
surrogate variable analysis. PLoS Genetics, 3(9):e161, 2007.

[11] Eugene Andres Houseman, John Molitor, and Carmen J Marsit. Reference-free cell mixture
adjustments in analysis of DNA methylation data. Bioinformatics, 30(10):1431-1439,
2014.

[12] Xiangyu Luo, Can Yang, and Yingying Wei. Detection of Cell-Type-Specific CpG Sites in Epigenome-Wide Association Studies.
