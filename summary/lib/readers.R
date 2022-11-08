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


read.invitroResults <- function() {
  dfInvitro = read.dbTables("invitro_selection_results")
  # dfInvitro = dfInvitro %>%
  #   mutate(
  #     ab_name=rx_name
  #   ) %>%
  #   inner_join(
  #     read.antibodies(),
  #     by = c("ab_name"),
  #     suffix = c(".invitro", ".mab")
  #   )
  dfInvitro
}

read.invivoResults <- function() {
  dfInvivo = read.dbTables("invivo_selection_results")
  dfInvivo
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

  dfSusc = dfSusc %>%
    mutate(
      ineffective = ifelse(is.na(ineffective), "", ineffective)
    ) %>%
    mutate(
      fold_cmp = ifelse(ineffective == 'experimental', ">", fold_cmp),
      fold = ifelse(ineffective == 'experimental', 100, fold)
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
        "resistant"
      )
    )
  )
  dfSusc$resistance_level = factor(
    dfSusc$resistance_level,
    levels = c("susceptible", "partial-resistance", "resistant")
  )
  dfSusc
}

read.suscResultsPlasam <- function(
  partialResistFold = 3,
  resistFold = 10,
  withExcluded = FALSE
) {
  dfSusc = read.suscResults(
    partialResistFold = partialResistFold,
    resistFold = resistFold,
    withExcluded = withExcluded
  )

  dfPlasma = read.dbTables(
      "rx_plasma",
      colClasses = c(
        titer = "character",
        subject_name = "character"
      )
  )
  if (withExcluded) {
    dfPlasma = bind_rows(dfPlasma, read.dbTables(
      "rx_plasma",
      colClasses = c(
        titer = "character",
        subject_name = "character"
      ),
      tables_dir = EXCLUDES_DIR
    ))
  }

  dfSubjectHistory = read.subjectHistory(withExcluded) %>%
    rename(
      infection=iso_name
    )

  dfPlasma = merge(
      dfPlasma,
      dfSubjectHistory,
      by = c("ref_name", 'subject_name')
    )

  dfCPIsolation = dfPlasma %>%
    filter(
      event == 'isolation',
      event_date == collection_date
    )

  dfCPInfection = dfPlasma %>%
    filter(
      event == 'infection'
    )

  dfCP = dfCPIsolation %>%
    left_join(
      dfCPInfection,
      by = c("ref_name", 'subject_name')
    ) %>%
    filter(
      (event.x == 'isolation' & event.y == 'infection')
      | (event.x == 'isolation' & is.na(event.y) )
    ) %>%
    select(
      ref_name,
      rx_name = rx_name.x,
      cumulative_group = cumulative_group.x,
      infection = infection.x
    ) %>%
    mutate(
      vaccine_name = ""
    ) %>%
    distinct()

  dfPlasma = merge(
      dfPlasma,
      dfSubjectHistory,
      by = c("ref_name", 'subject_name')
    )

  dfVP_1st = dfPlasma %>%
    filter(
      event.x == '1st dose isolation',
      event.y == '1st dose',
      event_date.x == collection_date
    )

  dfVP_2nd = dfPlasma %>%
    filter(
      event.x == '2nd dose isolation',
      event.y == '2nd dose',
      event_date.x == collection_date
    )
  dfVP_3rd = dfPlasma %>%
    filter(
      event.x == '3rd dose isolation',
      event.y == '3rd dose',
      event_date.x == collection_date
    )

  dfVP = bind_rows(dfVP_1st, dfVP_2nd, dfVP_3rd) %>%
    select(
      ref_name,
      rx_name,
      cumulative_group,
      vaccine_name = vaccine_name.y,
    ) %>%
    mutate(
      infection = ""
    ) %>%
    distinct()

  dfPlasma = bind_rows(dfCP, dfVP)

  merge(dfSusc, dfPlasma, by = c("ref_name", "rx_name"))
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

read.subjectHistory <- function(withExcluded) {
  dfSubjectHistory = read.dbTables(
      "subject_history",
      colClasses = c(
        titer = "character",
        subject_name = "character"
        )
      )

  if (withExcluded) {
    dfSubjectHistory = bind_rows(dfSubjectHistory, read.dbTables(
      "subject_history",
      colClasses = c(
        titer = "character",
        subject_name = "character"
      ),
      tables_dir = EXCLUDES_DIR
    ))
  }

  dfSubjectHistory
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
    mutate(
      ab_author_target = paste(target, collapse = ";")
    )

  dfMabPDBTargets = dfMAbTargets %>%
    filter(
      source == "structure"
    ) %>%
    group_by(ab_name) %>%
    mutate(
      ab_author_target = paste(target, collapse = ";")
    )

  dfMabTargets = bind_rows(dfMAbAuthorTargets, dfMabPDBTargets)

  read.dbTable("antibodies.csv") %>%
    left_join(dfMabTargets, by = "ab_name") %>%
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

read.variantMutations = function() {
  read.dbTable("isolate_mutations.csv")
}

read.virusVariants <- function() {
  dfMuts = read.variantMutations()
  dfSpikeMuts = dfMuts %>%
    filter(gene == "S") %>%
    mutate(mutation = paste0(position, amino_acid, split="")) %>%
    group_by(iso_name) %>%
    mutate(spike_muts=paste0(mutation, collapse="+")) %>%
    select(iso_name, gene, spike_muts) %>%
    unique()
  dfVariants = read.dbTable("isolates.csv")
  dfVariants = dfVariants %>%
    left_join(
      aggregate(
        cbind(num_muts = position) ~ iso_name,
        dfMuts,
        FUN = length
      ),
    by = c("iso_name"),
    suffix = c(".variant", ".mut")
    ) %>%
    mutate(
      num_muts = ifelse(is.na(num_muts), 0, num_muts)
    )
  dfVariants = dfVariants %>%
    left_join(
      dfSpikeMuts,
      by = c("iso_name"),
      suffix = c(".variant", ".mut")
    ) %>%
    mutate(
      spike_muts = ifelse(is.na(spike_muts), "", spike_muts)
    )
}
