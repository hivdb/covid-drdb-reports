library(DT)
library(slider)
library(dplyr)
library(kableExtra)

styledDT <- function(data, pageLength = 20) {
  data %>% datatable(
    escape = FALSE,
    rownames = FALSE,
    height = "auto",
    extensions = 'Buttons',
    options = list(pageLength = pageLength, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel', 'print'))
  )
}
styledKBL <- function(data, caption = NA) {
  options(knitr.kable.NA = '')
  data %>% kbl(caption = caption, row.names = FALSE) %>% kable_minimal(full_width = FALSE)
}

simpleOrdinal <- function(x) {
  # since scales::ordinal doesn't handle decimals
  x = as.character(x)
  if (endsWith(x, "1")) {
    x = sprintf("%sst", x)
  }
  else if (endsWith(x, "2")) {
    x = sprintf("%snd", x)
  }
  else if (endsWith(x, "3")) {
    x = sprintf("%srd", x)
  }
  else {
    x = sprintf("%sth", x)
  }
  return(x)
}

datatable.coverageDensity <- function(data, cov_breaks = c(0, 100, 300, 1000), displayPercentiles = c("0", "0.5", "1", "2", "5", "10")) {
  data = data[data$percentile %in% displayPercentiles,]
  for (cov_range in slide(cov_breaks, ~.x, .after = 1, .complete = TRUE)) {
    if (is.null(cov_range)) {
      break
    }
    data[between(data$coverage, cov_range[1], cov_range[2]), "coverageBin"] = sprintf(
      "%d-%d", cov_range[1], cov_range[2]
    )
    
  }
  data = aggregate(cbind("# samples" = sample_name) ~ coverageBin + percentile, data, FUN = length)
  dfList = lapply(displayPercentiles, function(x) {
    df = data[data$percentile == x, c("coverageBin", "# samples")]
    colnames(df) = c("Coverage", simpleOrdinal(x))
    covLevels = sort(unique(df$Coverage))
    df$Coverage = factor(df$Coverage, levels=covLevels)
    df
  })
  Reduce(
    function(x, y) merge(x = x, y = y, by = "Coverage", all = TRUE),
    dfList
  )
}

datatable.mutationIQR <- function(data, col_expr) {
  col_expr = eval(substitute(col_expr), data, parent.frame())
  data = Reduce(
    function(x, y) merge(x, y, by = "cutoff", all = TRUE),
    list(
      aggregate(cbind(median = col_expr) ~ cutoff, data, FUN=median),
      aggregate(cbind(p25 = col_expr) ~ cutoff, data, FUN=function(x) quantile(x, 0.25)),
      aggregate(cbind(p75 = col_expr) ~ cutoff, data, FUN=function(x) quantile(x, 0.75))
    )
  )
  data[order(data$cutoff),]
}

datatable.lowCovVariants <- function(dfVariants) {
  df = dfVariants[dfVariants$coverage < 1000,]
  df = df[df$amino_acid != "X",]  # remove all out-frame deletion
  df = df[order(df$is_vela_target, df$prevalence, decreasing =TRUE), c(
    "sample_name",
    "gene",
    "position",
    "ref_amino_acid",
    "amino_acid",
    "is_vela_target",
    "occurrence",
    "coverage",
    "prevalence",
    "drug_class"
  )]
  colnames(df) = c(
    "Sample", "Gene", "Pos", "RefAA", "VarAA", "VelaTarget",
    "VarCov", "PosCov", "Pcnt","DRM")
  df
}

datatable.discordDRMs <- function(dfVariantDRMs, cutoff) {
  df = dfVariantDRMs[dfVariantDRMs$prevalence >= cutoff | dfVariantDRMs$is_vela_target == TRUE,]
  df = df[df$prevalence < cutoff | df$is_vela_target != TRUE,]
  df$HIVDBCalled = ifelse(df$prevalence > cutoff, "Yes", "")
  df$VelaCalled = ifelse(df$is_vela_target == TRUE, "Yes", "")
  df = df[order(df$sample_name, -df$prevalence),
          c("sample_name", "gene", "position", "ref_amino_acid", "amino_acid",
            "HIVDBCalled", "VelaCalled", "occurrence", "prevalence", "drug_class")]
  colnames(df) = c("Sample", "Gene", "Pos", "RefAA", "VarAA", "HIVDBCalled", "VelaCalled", "VarCov", "Pcnt", "DRM")
  df$Gene = factor(df$Gene, levels = c("PR", "RT", "IN"))
  df[order(df$Sample, df$Gene, df$Pos, df$VarAA),]
}

datatable.intpnLevelDiffs <- function(data, dfDiscordDRMs, dfCutoffStats, cutoff) {
  cutoffText = sprintf("%g%%", cutoff * 100)
  dfDRMForPlot = dfDiscordDRMs %>%
    mutate(Mutation = ifelse(
      Pcnt < cutoff,
      sprintf("%s:%s%d%s†", Gene, RefAA, Pos, VarAA),
      sprintf("%s:%s%d%s", Gene, RefAA, Pos, VarAA)
    ))
  dfDRMForPlot = aggregate(
    cbind(HIVDBMuts = Mutation) ~ Sample, dfDRMForPlot[dfDRMForPlot$HIVDBCalled == "Yes",],
    FUN = function(x) paste0("+", x, collapse=", ")
  ) %>% merge(
    aggregate(
      cbind(VelaMuts = Mutation) ~ Sample, dfDRMForPlot[dfDRMForPlot$VelaCalled == "Yes",],
      FUN = function(x) paste0("-", x, collapse=", ")
    ),
    by = "Sample",
    all = TRUE
  )
  
  dfPlot = data %>%
    select(Sample, Drug, VelaLevel, HIVDBLevel) %>%
    mutate(
      VelaLevel = as.factor(substr(VelaLevel, 7, 7)),
      HIVDBLevel = as.factor(substr(HIVDBLevel, 7, 7))
    ) %>%
    mutate(
      "Vela => HIVDB" = as.integer(HIVDBLevel) - as.integer(VelaLevel)
    ) %>%
    merge(dfDRMForPlot, by = "Sample", all.x = TRUE) %>%
    merge(dfCutoffStats[dfCutoffStats$cutoff == cutoffText,], by.x = "Sample", by.y = "sample_name")

  dfAbsDiff = aggregate(
    cbind(maxAbsDiff = abs(`Vela => HIVDB`)) ~ Sample, dfPlot, FUN = max
  ) %>% merge(
    aggregate(
      cbind(meanAbsDiff = abs(`Vela => HIVDB`)) ~ Sample, dfPlot, FUN = mean
    ),
    by = "Sample"
  ) %>% merge(
    aggregate(
      cbind(minAbsDiff = abs(`Vela => HIVDB`)) ~ Sample, dfPlot, FUN = min
    ),
    by = "Sample"
  ) %>% merge(
    aggregate(
      cbind(medianAbsDiff = abs(`Vela => HIVDB`)) ~ Sample, dfPlot, FUN = median
    ),
    by = "Sample"
  )
  dfAbsDiff = dfAbsDiff[order(dfAbsDiff$meanAbsDiff),]
  dfAbsDiff$evenRow = seq(1, nrow(dfAbsDiff)) %% 2 == 0

  dfPlot = merge(dfPlot, dfAbsDiff, by = "Sample")
  dfPlot$Sample = factor(dfPlot$Sample, levels = dfAbsDiff$Sample)
  dfPlot
}