library(dplyr)
library(scales)

`%notin%` <- Negate(`%in%`)

.getpwd <- function() {
  pwd = getSrcDirectory(.getpwd)[1]
  if (is.na(pwd) || startsWith(pwd, '.')) {
    pwd = '.'
  }
  print(pwd)
  return(pwd)
}

TABLES_DIR = normalizePath(sprintf('%s/../covid-drdb-payload/tables', .getpwd()))
EXCLUDES_DIR = normalizePath(sprintf('%s/../covid-drdb-payload/excluded', .getpwd()))

read.dbTable <- function(filename, colClasses = NA, tables_dir = TABLES_DIR) {
  read.table(
    sprintf("%s/%s", tables_dir, filename),
    header = TRUE,
    sep = ",",
    quote = '"',
    na.strings = "NULL",
    colClasses = colClasses
  )
}

read.dbTables <- function(dirname, colClasses = NA, tables_dir = TABLES_DIR) {
  files = list.files(sprintf("%s/%s", tables_dir, dirname),
                     pattern = "*.csv", full.names = FALSE)
  if (length(files) == 0) {
    return(data.frame())
  }
  dfs = lapply(files, function(x) {
    read.dbTable(
      sprintf("%s/%s", dirname, x),
      colClasses = colClasses,
      tables_dir = tables_dir
    )
  })
  do.call(bind_rows, dfs)
}

read.suscResults <- function(
  partialResistFold = 3,
  resistFold = 10,
  withExcluded = FALSE
) {
  dfSusc = read.dbTables("susc_results", colClasses = c(fold = "numeric"))
  if (withExcluded) {
    dfSusc = dfSusc %>%
      mutate(include = TRUE) %>%
      bind_rows(
        read.dbTables(
          "susc_results",
          colClasses = c(fold = "numeric"),
          tables_dir = EXCLUDES_DIR
        ) %>%
        mutate(include = FALSE)
      )
  }
  dfSusc$fold_cmp = ifelse(
    !is.na(dfSusc$fold) &
    !is.null(dfSusc$fold) & (
      dfSusc$fold_cmp == "" |
      is.na(dfSusc$fold_cmp)
    ),
    yes = "=",
    no = dfSusc$fold_cmp)
  dfSusc$ordinal_number = ifelse(
    is.na(dfSusc$ordinal_number),
    yes = 1,
    no = dfSusc$ordinal_number
  )
  dfSusc$cumulative_count = ifelse(
    is.na(dfSusc$cumulative_count),
    yes = 1,
    no = dfSusc$cumulative_count
  )
  dfSusc$resistance_level = ifelse(
    !is.na(dfSusc$resistance_level),
    dfSusc$resistance_level,
    ifelse(
      dfSusc$fold_cmp %in% c("<", "~", "=") &
      dfSusc$fold <= partialResistFold,
      "susceptible",
      ifelse(
        dfSusc$fold_cmp %in% c("<", "~", "=") &
        dfSusc$fold <= resistFold,
        "partial-resistance",
        "resistance"
      )
    )
  )
  dfSusc$resistance_level = factor(
    dfSusc$resistance_level,
    levels = c("susceptible", "partial-resistance", "resistance")
  )
  dfSusc
}

read.suscResultsCP <- function(
  partialResistFold = 3,
  resistFold = 10,
  withExcluded = FALSE
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold,
    withExcluded = withExcluded
  )
  dfCP = read.dbTables("rx_conv_plasma")
  if (withExcluded) {
    dfCP = bind_rows(dfCP, read.dbTables(
      "rx_conv_plasma",
      tables_dir = EXCLUDES_DIR
    ))
  }
  merge(dfSusc, dfCP, by = c("ref_name", "rx_name"))
}

read.suscResultsIP <- function(
  partialResistFold = 3,
  resistFold = 10,
  withExcluded = FALSE
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold,
    withExcluded = withExcluded
  )
  dfIP = read.dbTables("rx_immu_plasma")
  if (withExcluded) {
    dfIP = bind_rows(dfIP, read.dbTables(
      "rx_immu_plasma",
      tables_dir = EXCLUDES_DIR
    ))
  }
  merge(dfSusc, dfIP, by = c("ref_name", "rx_name"))
}

read.suscResultsMAb <- function(
  partialResistFold = 3,
  resistFold = 10,
  withExcluded = FALSE
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold,
    withExcluded = withExcluded
  )
  dfMAb = read.dbTables("rx_antibodies")
  if (withExcluded) {
    dfMAb = bind_rows(dfMAb, read.dbTables(
      "rx_antibodies",
      tables_dir = EXCLUDES_DIR
    ))
  }
  dfMAb = dfMAb %>%
    inner_join(
      read.antibodies(),
      by = "ab_name",
      suffix = c(".rxmab", ".mab")
    )
  dfMAb = aggregate(cbind(
    ab_name,
    short_ab_name,
    short_ab_name_with_author_target
  ) ~ ref_name + rx_name, dfMAb, FUN = function(x) x)
  dfSusc %>%
    inner_join(
      dfMAb,
      by = c("ref_name", "rx_name"),
      suffix = c(".susc", ".rxmab")
    )
}

read.antibodies <- function() {
  dfMAbTargets = read.dbTable("antibody_targets.csv")
  dfMAbAuthorTargets = dfMAbTargets %>%
    filter(
      source == "author" &
      ab_name %notin% filter(
        dfMAbTargets,
        source == 'structure'
      )$ab_name
    ) %>%
    group_by(ab_name) %>%
    summarise(ab_author_target = paste(target, collapse = ";"))
  
  read.dbTable("antibodies.csv") %>%
    left_join(dfMAbAuthorTargets, by = "ab_name") %>%
    mutate(
      short_ab_name = ifelse(
        is.na(abbreviation_name),
        ab_name,
        abbreviation_name
      ),
      short_ab_name_with_author_target = ifelse(
        is.na(ab_author_target),
        short_ab_name,
        sprintf("%s (%s)", short_ab_name, ab_author_target)
      )
    )
}

read.articles <- function() {
  read.dbTable("articles.csv")
}

read.strainMutations = function() {
  read.dbTable("strain_mutations.csv")
}

read.virusStrains <- function() {
  dfMuts = read.strainMutations()
  dfStrains = read.dbTable("virus_strains.csv")
  merge(
    dfStrains,
    aggregate(
      cbind(num_muts = position) ~ strain_name,
      dfMuts,
      FUN = length
    ),
    by = "strain_name"
  )
}