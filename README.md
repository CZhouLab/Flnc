# Flnc

## Introduction

Flnc is software that can accurately identify full-length long noncoding RNAs (lncRNAs) from human RNA-seq data. lncRNAs are linear transcripts of more than 200 nucleotides that do not encode proteins. The most common approach for identifying lncRNAs from RNA-seq data which examines the coding abilities of assembled transcripts will result in a very high false-positive rate (30%-75%) of lncRNA identification. The falsely discovered lncRNAs lack transcriptional start sites and most of them are RNA fragments or result from transcriptional noise. Unlike the false-positive lncRNAs, true lncRNAs are full-length lncRNA transcripts that include transcriptional start sites (TSSs). To exclude these false lncRNAs, H3K4me3 chromatin immunoprecipitation sequencing (ChIP-seq) data had been used to examine transcriptional start sites of putative lncRNAs, which are transcripts without coding abilities. However, because of cost, time, and the limited availability of sample materials for generating H3K4me3 ChIP-seq data, most samples (especially clinical biospecimens) may have available RNA-seq data but lack matched H3K4me3 ChIP-seq data. This Flnc method solves the problem of lacking transcriptional initiation profiles when identifying lncRNAs. 

Flnc integrates seven machine-learning algorithms built with four genomic features. Flnc achieves state-of-the-art prediction power with a AUROC score over 0.92. Flnc significantly improves the prediction accuracy from less than 50% using the common approach to over 85% on five independent datasets without requiring matched H3K4me3 ChIP-seq data. In addition to the stranded polyA-selected RNA-seq data, Flnc can also be applied to identify lncRNAs from stranded RNA-seq data of ribosomal RNA depleted samples or unstranded RNA-seq data of polyA-selected samples. 

![workflow](Picture1.png)

**Please cite our paper at BioRxiv, if you find Flnc useful for your research. The paper has been submitted to a peer-reviewed journal**

Version: 1.0.0

Last Modified: 07/27/2022

Authors: Zixiu Li (zixiu.li@umassmed.edu), Chan Zhou (chan.zhou@umassmed.edu)

Maintainer: Zixiu Li


## Prerequisites

To use Flnc, you will need the following programs in your PATH:

•       singularity (>=3.7.1)

•       python2 (>=2.7.18) 

•       OS: high performance computing cluster in Linux (suggested)

•       Reference genome: hg38

