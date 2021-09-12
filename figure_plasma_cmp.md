---
title: "Susceptibility summary figures"
output:
  html_document:
    keep_md: yes
    css: "css/main.css"
---

Report generated at: 09/12/2021 11:47 AM PDT






## Plasma comparison {.tabset}

### Fold
<img src="figure_plasma_cmp_files/figure-html/fold-cmp-1m-1.png" width="1920" />



### Titer
<img src="figure_plasma_cmp_files/figure-html/titer-cmp-1m-1.png" width="1920" />



### Fold vs Titer
<img src="figure_plasma_cmp_files/figure-html/fold-vs-titer-1m-1.png" width="1920" />


<!-- ### Titer 1M vs Titer other -->
<!-- ```{r titer-1m-vs-other, fig.height=20, fig.retina=4, fig.width=20, message=FALSE, warning=FALSE, echo=FALSE} -->

<!-- titer_df = read.dbTable('figure_plasma_titer.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='1') -->
<!-- titer_other = read.dbTable('figure_plasma_titer.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='>=2') -->

<!-- plots = list() -->
<!-- i = 1 -->
<!-- line_width = 4 -->
<!-- j = i + line_width -->
<!-- rx_name_length = length(RX_NAME_LIST) -->

<!-- for (rx in RX_NAME_LIST) { -->
<!--   dfPlot = titer_df %>% -->
<!--     filter(rx_name == rx) %>% -->
<!--     filter(variant %in% VARIANT) %>% -->
<!--     mutate( -->
<!--       variant = factor(variant, VARIANT), -->
<!--       level = factor(level, TITER_LEVEL) -->
<!--     ) -->

<!--   if (rx == 'CP') { -->
<!--     figLabel = "Convalescent plasma" -->
<!--     dfPlot = dfPlot %>% -->
<!--       filter(infection == 'infected') -->
<!--   } -->
<!--   else { -->
<!--     figLabel = rx -->
<!--     dfPlot = dfPlot %>% -->
<!--       filter(infection == 'naive') -->
<!--   } -->

<!--   p = draw_rx_titer(dfPlot, figLabel) -->

<!--   if (i %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--     ) + -->
<!--       ylab("# Results (1M)") -->
<!--   } -->

<!--   plots[[i]] = p -->
<!--   i = i + 1 -->

<!--   dfPlot = titer_other %>% -->
<!--     filter(rx_name == rx) %>% -->
<!--     filter(variant %in% VARIANT) %>% -->
<!--     mutate( -->
<!--       variant = factor(variant, VARIANT), -->
<!--       level = factor(level, TITER_LEVEL) -->
<!--     ) -->
<!--   p = draw_rx_titer(dfPlot, figLabel) -->

<!--   if (j %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--     ) + -->
<!--       ylab("# Results (>=2M)") -->
<!--   } -->

<!--   if (j > (2 * rx_name_length - line_width)) { -->
<!--     p = p + theme( -->
<!--       axis.title.x = element_blank(), -->
<!--       axis.text.x = element_text(size=14, angle=45, hjust=1), -->
<!--     ) -->
<!--   } -->

<!--   plots[[j]] = p -->
<!--   j = i + line_width -->
<!-- } -->

<!-- grid.arrange(grobs = plots, ncol=line_width) -->
<!-- # grid.arrange(grobs = plots) -->
<!-- ``` -->


<!-- ### Fold infected vs uninfected -->
<!-- ```{r fold-infect-cmp, fig.height=20, fig.retina=4, fig.width=20, message=FALSE, warning=FALSE, echo=FALSE} -->

<!-- fold_df = read.dbTable('figure_plasma_fold.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='1') -->

<!-- plots = list() -->
<!-- i = 1 -->
<!-- line_width = 4 -->
<!-- j = i + line_width -->
<!-- rx_name_length = length(RX_NAME_LIST) - 1 -->

<!-- for (rx in RX_NAME_LIST) { -->
<!--     if (rx == 'CP') { -->
<!--       next -->
<!--     } -->

<!--     dfPlot = fold_df %>% -->
<!--       filter(rx_name == rx) %>% -->
<!--       filter(variant %in% VARIANT) %>% -->
<!--       mutate( -->
<!--         variant = factor(variant, VARIANT), -->
<!--         susc = factor(susc, SUSC_VALUE) -->
<!--       ) -->
<!--     figLabel = rx -->

<!--     dfPlotNaive = dfPlot %>% -->
<!--       filter(infection == 'naive') -->

<!--     p = draw_rx_fold(dfPlotNaive, figLabel) -->

<!--     if (i %% line_width == 1) { -->
<!--       p = p + theme( -->
<!--         axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--         ) + -->
<!--       ylab("# Results (naive)") -->
<!--     } -->

<!--     plots[[i]] = p -->
<!--     i = i+ 1 -->

<!--     dfPlotInfected = dfPlot %>% -->
<!--       filter(infection == 'infected') -->

<!--     p = draw_rx_fold(dfPlotInfected, figLabel) -->

<!--     if (j %% line_width == 1) { -->
<!--       p = p + theme( -->
<!--         axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--       ) + -->
<!--         ylab("# Results (infected)") -->
<!--     } -->

<!--     if (j > (2 * rx_name_length - line_width)) { -->
<!--       p = p + theme( -->
<!--         axis.title.x = element_blank(), -->
<!--         axis.text.x = element_text(size=14, angle=45, hjust=1), -->
<!--       ) -->
<!--     } -->

<!--     plots[[j]] = p -->
<!--     j = i + line_width -->
<!-- } -->

<!-- while (i %% 4 != 1) { -->
<!--   p = draw_blank() -->
<!--   plots[[i]] = p -->
<!--   i = i + 1 -->
<!-- } -->

<!-- grid.arrange(grobs = plots, ncol=line_width) -->
<!-- # grid.arrange(grobs = plots) -->
<!-- ``` -->


<!-- ### Titer infected vs uninfected -->
<!-- ```{r titer-infect-cmp, fig.height=20, fig.retina=4, fig.width=20, message=FALSE, warning=FALSE, echo=FALSE} -->

<!-- titer_df = read.dbTable('figure_plasma_titer.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='1') -->

<!-- plots = list() -->
<!-- i = 1 -->
<!-- line_width = 4 -->
<!-- j = i + line_width -->
<!-- rx_name_length = length(RX_NAME_LIST) - 1 -->

<!-- for (rx in RX_NAME_LIST) { -->
<!--   if (rx == 'CP') { -->
<!--     next -->
<!--   } -->

<!--   dfPlot = titer_df %>% -->
<!--     filter(rx_name == rx) %>% -->
<!--     filter(variant %in% VARIANT) %>% -->
<!--     mutate( -->
<!--       variant = factor(variant, VARIANT), -->
<!--       level = factor(level, TITER_LEVEL) -->
<!--     ) -->
<!--   figLabel = rx -->

<!--   dfPlotNaive = dfPlot %>% -->
<!--     filter(infection == 'naive') -->

<!--   p = draw_rx_titer(dfPlotNaive, figLabel) -->

<!--   if (i %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--     ) + -->
<!--       ylab("# Results (naive)") -->
<!--   } -->

<!--   plots[[i]] = p -->
<!--   if ((i + 1) %% 8 <= 4) { -->
<!--     i = i + 1 -->
<!--   } else { -->
<!--     i = i + 1 + line_width -->
<!--   } -->

<!--   dfPlotInfected = dfPlot %>% -->
<!--     filter(infection == 'infected') -->

<!--   p = draw_rx_titer(dfPlotInfected, figLabel) -->

<!--   if (j %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--     ) + -->
<!--       ylab("# Results (infected)") -->
<!--   } -->

<!--   if (j > (2 * rx_name_length - line_width)) { -->
<!--     p = p + theme( -->
<!--       axis.title.x = element_blank(), -->
<!--       axis.text.x = element_text(size=14, angle=45, hjust=1), -->
<!--     ) -->
<!--   } -->

<!--   plots[[j]] = p -->
<!--   j = i + line_width -->
<!-- } -->

<!-- while (i %% 4 != 1) { -->
<!--   p = draw_blank() -->
<!--   plots[[i]] = p -->
<!--   i = i + 1 -->
<!-- } -->

<!-- grid.arrange(grobs = plots, ncol=line_width) -->
<!-- # grid.arrange(grobs = plots) -->
<!-- ``` -->

<!-- ## Titer 1M violin plot -->
<!-- ```{r titer-1m-violin, fig.height=20, fig.retina=4, fig.width=20, message=FALSE, warning=FALSE, echo=FALSE} -->

<!-- titer_df = read.dbTable('figure_plasma_titer_points.csv', tables_dir=TABLEDIR) %>% -->
<!--   filter(month=='1') -->

<!-- plots = list() -->
<!-- i = 1 -->
<!-- line_width = 4 -->
<!-- j = i + line_width -->
<!-- rx_name_length = length(RX_NAME_LIST) -->

<!-- for (rx in RX_NAME_LIST) { -->
<!--   dfPlot = titer_df %>% -->
<!--     filter(rx_name == rx) %>% -->
<!--     filter(variant %in% VARIANT) %>% -->
<!--     mutate( -->
<!--       variant = factor(variant, VARIANT) -->
<!--     ) -->

<!--   if (rx == 'CP') { -->
<!--     figLabel = "Convalescent plasma" -->
<!--     dfPlot = dfPlot %>% -->
<!--       filter(infection == 'infected') -->
<!--   } -->
<!--   else { -->
<!--     figLabel = rx -->
<!--     dfPlot = dfPlot %>% -->
<!--       filter(infection == 'naive') -->
<!--   } -->

<!--   p = draw_titer_violin(dfPlot, figLabel) -->

<!--   if (i %% line_width == 1) { -->
<!--     p = p + theme( -->
<!--       axis.title.y = element_text(size=15, face="bold", colour = "black"), -->
<!--     ) + -->
<!--       ylab("NT50") -->
<!--   } -->

<!--   if (i > rx_name_length - line_width) { -->
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



## Titer 1M compare by variant
<img src="figure_plasma_cmp_files/figure-html/titer-1m-variant-cmp-1.png" width="1920" />


## Titer 1M compare by Plasma
<img src="figure_plasma_cmp_files/figure-html/titer-1m-cmp-plasma-boxplot-1.png" width="1920" />


## Titer 1M vs 6M boxplot
<img src="figure_plasma_cmp_files/figure-html/titer-1m-vs-6m-boxplot-1.png" width="1920" />


## Titer VP Naive vs Infection boxplot
<img src="figure_plasma_cmp_files/figure-html/titer-naive-vs-infection-boxplot-1.png" width="1920" />
