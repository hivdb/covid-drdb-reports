---
title: "Omicron Mab Tree figure"
output:
  html_document:
    keep_md: yes
---

Report generated at: 01/31/2022 05:28 PM PDT






## Figures {.tabset}

### BAM and ETE

```
## [1] "BAM"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 48, p-value = 0.05731
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "BAM"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 30, p-value = NA
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "ETE"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 36, p-value = 0.2065
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "ETE"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 20, p-value = 0.2031
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "BAM/ETE"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 10, p-value = 0.2286
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "BAM/ETE"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 6, p-value = NA
## alternative hypothesis: true location shift is not equal to 0
```

<img src="Omicron_mAb_compare_virus_type_files/figure-html/bam-ete-1.png" width="768" />

### CAS and IMD

```
## [1] "CAS"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 44.5, p-value = 0.2686
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "CAS"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 31, p-value = 0.7878
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "IMD"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 48.5, p-value = 0.1312
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "IMD"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 33, p-value = NA
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "CAS/IMD"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 11.5, p-value = 0.8057
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "CAS/IMD"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 14, p-value = 0.2404
## alternative hypothesis: true location shift is not equal to 0
```

<img src="Omicron_mAb_compare_virus_type_files/figure-html/cas-imd-1.png" width="768" />

### SOT

```
## [1] "SOT"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 57, p-value = 0.1042
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "SOT"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 63, p-value = 0.02677
## alternative hypothesis: true location shift is not equal to 0
```

<img src="Omicron_mAb_compare_virus_type_files/figure-html/sot-1.png" width="768" />

### CIL and TIX

```
## [1] "CIL"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 53, p-value = 0.002539
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "CIL"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 11.5, p-value = 0.07446
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "TIX"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 51, p-value = 0.005572
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "TIX"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 28.5, p-value = 0.9053
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "CIL/TIX"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 19, p-value = 0.03175
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "CIL/TIX"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 11, p-value = 0.9048
## alternative hypothesis: true location shift is not equal to 0
```

<img src="Omicron_mAb_compare_virus_type_files/figure-html/cil-tix-1.png" width="768" />

### Other

```
## [1] "ADI"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 2, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "ADI"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum exact test
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 2, p-value = 1
## alternative hypothesis: true location shift is not equal to 0
```

```
## [1] "REG"
## [1] "Wildtype: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_wt and pseudovirus_wt
## W = 10, p-value = 0.2118
## alternative hypothesis: true location shift is not equal to 0
## 
## [1] "REG"
## [1] "Omicron: pseudovirus IC50 < live virus IC50"
## 
## 	Wilcoxon rank sum test with continuity correction
## 
## data:  live_virus_omicron and pseudovirus_omicron
## W = 6, p-value = NA
## alternative hypothesis: true location shift is not equal to 0
```

```
## NULL
```

```
## NULL
```

```
## NULL
```

```
## NULL
```

```
## NULL
```

```
## NULL
```

<img src="Omicron_mAb_compare_virus_type_files/figure-html/other-1.png" width="768" />
