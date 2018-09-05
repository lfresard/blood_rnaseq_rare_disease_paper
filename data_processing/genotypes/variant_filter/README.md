# udn / variant filtering pipeline

The variant QC generally follows the ExAC guideline as in [ref1](https://doi.org/10.1016/j.ajhg.2018.05.002) and [ref2](https://doi.org/10.1038/nature19057).

### variant based metrics:
#### 1. keep only "PASS" variants, determined by the GATK variant quality score recalibration (VQSR)

adhoc filter for non-GATK processed vcfs: keep "." or "num prev seen samples > 30"

#### 2. as in ["hail"](https://hail.is/docs/stable/index.html) api

```
.filter_genotypes(
         """
         (g.isHomRef && (g.ad[0] / g.dp < 0.8 || g.gq < 20 || g.dp < 20)) ||
         (g.isHomVar && (g.ad[1] / g.dp < 0.8 || g.pl[0] < 20 || g.dp < 20)) ||
         (g.isHet && ( (g.ad[0] + g.ad[1]) / g.dp < 0.8 || g.ad[1] / g.dp < 0.20 || g.pl[0] < 20 || g.dp < 20))
         """, keep = False)
```         


as in VCF filters:

Filters | FORMAT/GT | FORMAT/DP | FORMAT/GQ | FORMAT/PL | FORMAT/AD
---|---|---|---|---|---
keep|=='0/0'|>20|>20|PL[0]<20|AD[0]/DP>0.8
keep|=='0/1'|>20|>20|PL[1]<20|AD[0]/DP>0.2 & AD[1]/DP>0.2 & (AD[0]+AD[1])/DP>0.2
keep|=='1/1'|>20|>20|PL[2]<20|AD[1]/DP>0.8

note 1: PL = *-10\*log10(likelihood)*, most likely being *PL=0*;  
note 2: FORMAT information is different across sits/samples/institutions, cannot apply uniform filters
* CGS	GT:AD:DP:GQ:PL
* CGS	GT:GQ:PL
* CHEO	GT
* CHEO	GT:AD:DP:GQ:PL
* CHEO	GT:GQ:PL
* CHEO	GT:PL
* CHEO	GT:PL:DP:SP:GQ
* UDN	GT:AD
* UDN	GT:AD:DP
* UDN	GT:AD:DP:GQ:PGT:PID:PL
* UDN	GT:AD:DP:GQ:PL
* UDN	GT:AD:DP:PGT:PID
* UDN	GT:AD:GQ:PGT:PID:PL
* UDN	GT:AD:GQ:PL
* UDN	GT:AD:PGT:PID
* UDN	GT:VR:RR:DP:GQ

### site based metrics: 
#### 3. exclude variants with HWE p-value < 1e-6
#### 4. exclude varriants with call rate < 0.80 (missing > 20%)

[ref1](https://doi.org/10.1016/j.ajhg.2018.05.002) Ganna, A., Satterstrom, F. K., Zekavat, S. M., Das, I., Kurki, M. I., Churchhouse, C., … Neale, B. M. (2018). Quantifying the Impact of Rare and Ultra-rare Coding Variation across the Phenotypic Spectrum. American Journal of Human Genetics, 102(6), 1204–1211.  https://github.com/andgan/ultra_rare_pheno_spectrum

[ref2](https://doi.org/10.1038/nature19057) Lek, M., Karczewski, K. J., Minikel, E. V, Samocha, K. E., Banks, E., Fennell, T., … Exome Aggregation Consortium. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature, 536(7616), 285–91. 
