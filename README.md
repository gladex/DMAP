### PanCanMAP: Pan-Cancer DNA Methylation Analysis Package for ENCODE

This GitHub portal is for depositing the Pan-Cancer DNA Methylation Analysis Package (PanCanMAP) & related analysis results for our recent ENCODE project.

#### D3_T47D.methyl.DVQV.25p:
Analysis and annotation results for the significant differentially-methylated CpG sites (SDMC) with reference to one cell type, here we selected T-47D as the study case. The results are further filtered based on the lifted methylation difference threshold (at least 25% methylation difference for the paired groups). And the SDMC list contains 106,252 DMCs, together the related statistical p-value and adjusted q-value are also provided; 

#### D4_T47D_SigDMR:
Statistical analysis and annotation results for the differentially-methylated regions (DMR) with reference to one cell type, for consistence we selected T-47D as the case. We identified 16,277 DMR candidates from all the DMCs, with the adjusted q-value <= 0.01, CpG base methylation difference cutoff, 25, and DMR mean methylation difference cutoff, 20. Within those candidates, 8,936 entries present hyper-methylated and 7,341 with hypo-methylated status. With the lifted thresholds, namely adjusted q-value <= 0.001, differentially-methylated CpG base count >= 5, we further detected 7,537 significant DMRs (Sig-DMRs), where 3,512 entries are significantly hypermethylated-DMRs (Sig-Hyper-DMRs), and 4,025 significantly hypomethylated-DMRs (Sig-Hypo-DMRs);

#### D5_T47D_SigDMR.HYPER:
Statistical analysis and annotation results for the significantly hypermethylated-DMRs (Sig-Hyper-DMRs) with reference to T-47D cell type;

#### D6_T47D_SigDMR.HYPO:
Statistical analysis and annotation results for the significantly hypomethylated-DMRs (Sig-Hypo-DMRs) with reference to T-47D cell type;

#### D7_T47D_methyGenes:
Statistical analysis and annotation results for the identified genes from all DMRs (hyper-DMRs and hypo-DMRs) with reference to T-47D cell type;

#### D8_GO_HyperDMR:
Gene Ontology annotation results for the TSGs identified from hyper-DMRs;

#### D9_GO_HypoDMR:
Gene Ontology annotation results for the TSGs identified from hypo-DMRs.

All the curated reference sources and read-to-use pan-cancer analysis results are deposited at this GitHub portal. We continuously update and revise the contents of this package and related results. Suggestions or advices to pan-cancer analysis repository, PanCanMAP, are very welcome to me at bh.tang@outlook.com.

Any citation or usage of the package, please cite: 

Tang, B.; Zhou, Y.; Wang, C.-M.; Huang, T.H.M.; Jin, V.X. Integration of DNA methylation and gene transcription across nineteen cell types reveals cell type-specific and genomic region-dependent regulatory patterns. Scientific Reports 2017, 7, 3626.

[EFG@HHU | BHT | Last update: 2017/04/15]
