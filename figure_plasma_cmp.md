---
title: "Susceptibility summary figures"
output:
  html_document:
    keep_md: yes
    css: "css/main.css"
---

Report generated at: 09/15/2021 05:48 PM PDT



## Number results {.tabset}

### Fold 1M
<img src="figure_plasma_cmp_files/figure-html/fold-cmp-1m-1.png" width="1920" />

### Titer 1M
<img src="figure_plasma_cmp_files/figure-html/titer-cmp-1m-1.png" width="1920" />


## Titer vs Fold {.tabset}

<!-- ### Titer percent -->
<!-- ```{r titer-percent, fig.height=20, fig.retina=4, fig.width=20, message=FALSE, warning=FALSE, echo=FALSE} -->

<!-- draw_titer_pcnt <- function(dfPlot, figLabel) { -->
<!--   dfAgg = dfPlot %>% -->
<!--     group_by(variant, num_study) %>% -->
<!--     summarise(num_results=sum(num_result)) %>% -->
<!--     mutate( -->
<!--       label = ifelse( -->
<!--         num_study > 1, -->
<!--         yes=sprintf('%s studies\n%s results', num_study, num_results), -->
<!--         no=ifelse( -->
<!--           num_study > 0, -->
<!--           yes=sprintf('%s study\n%s results', num_study, num_results), -->
<!--           no='' -->
<!--         ) -->
<!--       ) -->
<!--     ) -->

<!--   p <- ggplot(data=dfPlot, aes(x=variant, y=num_result)) + -->
<!--     geom_bar( -->
<!--       aes(fill=level), -->
<!--       position="fill", -->
<!--       stat="identity", -->
<!--       colour="black") + -->
<!--     scale_fill_manual( -->
<!--       values = c( -->
<!--         'low' = "#146aa8", -->
<!--         'middle' = "#7fcbee", -->
<!--         'high' = "#ffffff" -->
<!--       ), -->
<!--       name="NT titer", -->
<!--       breaks=TITER_LEVEL, -->
<!--       labels=c( -->
<!--         "<=40", -->
<!--         "40~100", -->
<!--         ">100" -->
<!--       ) -->
<!--     ) + -->
<!--     scale_y_continuous( -->
<!--       expand=expansion(mult = c(0, .1)), -->
<!--       breaks=c(0, 0.25, 0.5, 0.75, 1), -->
<!--       labels=c('0', '25%', '50%', '75%', '100%') -->
<!--     ) + -->
<!--     geom_text( -->
<!--       data=dfAgg, -->
<!--       aes(label=label), -->
<!--       size = 4, -->
<!--       vjust = -0.2, -->
<!--       y = 1 -->
<!--     ) + -->
<!--     theme( -->
<!--       legend.position = "none", -->
<!--       # legend.title = element_blank(), -->
<!--       axis.title.x = element_blank(), -->
<!--       axis.title.y = element_blank(), -->
<!--       axis.text.x = element_blank(), -->
<!--       axis.text.y = element_blank(), -->
<!--       plot.title = element_text(hjust = 0.5, size=12), -->
<!--       strip.background = element_blank(), -->
<!--       strip.text = element_blank(), -->
<!--       panel.background = element_blank(), -->
<!--       panel.grid.major = element_blank(), -->
<!--       panel.grid.minor = element_blank(), -->
<!--       panel.border = element_rect(fill = NA, color = 'black'), -->
<!--       plot.margin = margin(t=0, l=10, r=10, b=10) -->
<!--     ) -->

<!--   p = p + annotate( -->
<!--     "text", -->
<!--     y= Inf, -->
<!--     x = Inf, -->
<!--     label=figLabel, -->
<!--     size=5, -->
<!--     vjust=2, -->
<!--     hjust=1.1) -->
<!--   p -->
<!-- } -->

<!-- titer_df = read.dbTable('figure_plasma_titer.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='1') -->
<!-- titer_df = titer_df %>% -->
<!--   filter(variant %in% VARIANT) %>% -->
<!--   filter(rx_name %in% RX_NAME_LIST) -->
<!-- titer_df_CP = titer_df %>% -->
<!--   filter(rx_name == 'CP', infection == 'infected') -->
<!-- titer_df_VP = titer_df %>% -->
<!--   filter(rx_name != 'CP', infection =='naive') -->
<!-- titer_df = rbind(titer_df_VP, titer_df_CP) -->

<!-- plots = list() -->
<!-- i = 1 -->
<!-- line_width = 4 -->
<!-- num_graph = length(RX_NAME_LIST) - 1 -->

<!-- for (rx in RX_NAME_LIST) { -->
<!--   if (rx == 'CP') { -->
<!--     next -->
<!--   } -->
<!--   dfPlot = titer_df %>% -->
<!--     filter(rx_name == rx) %>% -->
<!--     mutate( -->
<!--       variant = factor(variant, VARIANT), -->
<!--       level = factor(level, TITER_LEVEL) -->
<!--     ) -->

<!--   figLabel = rx -->
<!--   p = draw_titer_pcnt(dfPlot, figLabel) -->

<!--   if (i %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--       axis.text.y = element_text(size=14) -->
<!--     ) + -->
<!--       ylab("% Results") -->
<!--   } -->

<!--   if (i > num_graph - line_width) { -->
<!--     p = p + theme( -->
<!--       axis.title.x = element_blank(), -->
<!--       axis.text.x = element_text(size=14, angle=45, hjust=1), -->
<!--     ) -->
<!--   } -->

<!--   plots[[i]] = p -->
<!--   i = i + 1 -->
<!-- } -->

<!-- grid.arrange(grobs = plots, ncol=line_width) -->

<!-- ``` -->

### Fold 1M vs Titer 1M
<img src="figure_plasma_cmp_files/figure-html/fold-vs-titer-1m-1.png" width="1920" />

### Fold 1M vs Titer 1M same scale
<img src="figure_plasma_cmp_files/figure-html/fold-vs-titer-1m-same-scale-1.png" width="1920" />


### Fold 1M vs Titer 1M Indiv
<img src="figure_plasma_cmp_files/figure-html/fold-vs-titer-1m-indiv-1.png" width="1920" />

## Titer by variant and by vaccine {.tabset}

### Titer 1M compare by variant
<img src="figure_plasma_cmp_files/figure-html/titer-1m-variant-cmp-1.png" width="1920" />

### Titer 1M compare by vaccine
<img src="figure_plasma_cmp_files/figure-html/titer-1m-cmp-plasma-boxplot-1.png" width="1920" />

## Titer fold correlation  {.tabset}

### Log axis
<img src="figure_plasma_cmp_files/figure-html/titer-fold-point-details-1.png" width="1920" />

### Log axis without one study
<img src="figure_plasma_cmp_files/figure-html/titer-fold-point-details-sp-1.png" width="1920" />

### Normal axis
<img src="figure_plasma_cmp_files/figure-html/titer-fold-point-details-normal-axis-1.png" width="1920" />

### Normal axis without one study
<img src="figure_plasma_cmp_files/figure-html/titer-fold-point-details-normal-axis-sp-1.png" width="1920" />

### Titer fold correlation by variant and VP
<img src="figure_plasma_cmp_files/figure-html/titer-fold-correlation-1.png" width="1920" />

## Other titer compare {.tabset}

### Titer 1M vs >=2M, NT50
<img src="figure_plasma_cmp_files/figure-html/titer-1m-vs-6m-boxplot-1.png" width="1920" />

### Titer 1M vs >=2M, num results
<img src="figure_plasma_cmp_files/figure-html/titer-1m-vs-other-1.png" width="1920" />


### Titer VP Naive vs Infection
<img src="figure_plasma_cmp_files/figure-html/titer-naive-vs-infection-boxplot-1.png" width="1920" />



### Titer VP Naive vs infection, num results
<img src="figure_plasma_cmp_files/figure-html/titer-infect-cmp-1.png" width="1920" />

### Titer 1M violin plot
<img src="figure_plasma_cmp_files/figure-html/titer-1m-violin-1.png" width="1920" />


## Other Fold compare

### Fold naive vs uninfected, num results
<img src="figure_plasma_cmp_files/figure-html/fold-infect-cmp-1.png" width="1920" />
