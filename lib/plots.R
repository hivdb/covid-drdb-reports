library(dplyr)
library(ggplot2)
library(gridExtra)

blank_theme <- function() {
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.position = "none",
    panel.border = element_rect(color = "black", size = 1, fill = "transparent")
  )
}

angled_x_axis <- function() {
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
}

plot.coverageDensity <- function(data, binwidth = 100, limits = c(-5, 1000), displayPercentiles = c("0", "0.5", "1", "2", "5", "10")) {
  major_tick = binwidth * 5
  maxBreak = ceiling(max(data$coverage) / major_tick) * major_tick
  p2 = ggplot(data[data$percentile %in% displayPercentiles,],
              aes(x = coverage, fill = percentile, order=desc(percentile))) +
    scale_y_continuous("# samples", expand = c(0, 0)) +
    blank_theme() + angled_x_axis()
  breaks = seq(from = 0, to = maxBreak, by = binwidth)
  # p1 = p2 +
  #   geom_histogram(alpha = .5, binwidth = binwidth, position = "identity") +
  #   geom_density(aes(y = ..count.. * binwidth), alpha = .4, bw = binwidth, position = "identity") +
  #   scale_x_continuous(
  #     "Read coverage",
  #     expand = c(0, 0),
  #     limits = limits,
  #     breaks = breaks,
  #     labels = ifelse(breaks %% major_tick == 0, breaks, "")
  #   ) +
  #   theme(axis.title.x = element_text(hjust = 2/3))
  p2 = p2 +
    geom_histogram(binwidth = binwidth) +
    scale_x_continuous(
      "Read coverage",
      expand = c(0, 1),
      breaks = breaks,
      labels = ifelse(breaks %% major_tick == 0, breaks, "")
    ) +
    coord_cartesian(xlim = limits) +
    facet_wrap(~ percentile, ncol = 1, strip.position = "right")
  # grid.arrange(p1, p2, layout_matrix = rbind(
  #  c(1, 1, 1, 2)
  # ))
  p2
}

plot.mutationDistribution <- function(data, col_expr, x_label, binwidth = 5) {
  col_expr = eval(substitute(col_expr), data, parent.frame())
  major_tick = binwidth * 5
  max_break = ceiling(max(col_expr) / major_tick) * major_tick
  breaks = seq(from = 0, to = max_break, by = binwidth)
  p2 = ggplot(data, aes(x = col_expr, fill = cutoff)) +
    blank_theme() + angled_x_axis()
  p1 = p2 +
    geom_histogram(position = "identity", binwidth = binwidth, alpha = 0.7) +
    # geom_density(aes(y = ..count.. * binwidth), position = "identity", bw = binwidth, alpha = 0.7) +
    scale_x_continuous(
      x_label,
      expand=c(0, 0),
      breaks = breaks,
      labels = ifelse(breaks %% major_tick == 0, breaks, "")
    ) +
    scale_y_continuous("# samples", expand=c(0, 0)) +
    theme(axis.title.x = element_text(hjust = 3/4))
  p2 = p2 +
    ggtitle("Cutoff") +
    geom_histogram(binwidth = binwidth) +
    scale_x_continuous("", expand=c(0, binwidth)) +
    scale_y_continuous(expand=c(0, 0)) +
    theme(
      axis.title.y = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5)
    ) +
    facet_wrap(~ cutoff, ncol = 1, strip.position = "right")
  grid.arrange(p1, p2, layout_matrix = rbind(
    c(1, 1, 1, 1, 1, 1, 2, 2, 2)
  ))
}

plot.intpnLevelDiffs <- function(data) {
  maxLenHIVDBMuts = max(nchar(data$HIVDBMuts), na.rm=TRUE)
  maxLenVelaMuts = max(nchar(data$VelaMuts), na.rm=TRUE)
  maxDiff = 4
  word2WidthScale = 0.2
  offsetNumUnusual = maxDiff + 1
  offsetNumApobec = offsetNumUnusual + 2.5
  offsetNumStop = offsetNumApobec + 2.5
  offsetVelaOnly = offsetNumStop + 2.5
  offsetHIVDBOnly = offsetVelaOnly + 1 + maxLenVelaMuts * word2WidthScale
  upperLimit = offsetHIVDBOnly + maxLenHIVDBMuts * word2WidthScale
  
  plot = ggplot(data, aes(x = Sample)) +
    ggtitle("Resistance Level Differences (-Vela, +HIVDB)") +
    geom_tile(
      mapping = aes(
        y = (upperLimit - maxDiff) / 2,
        width = 1,
        height = ifelse(evenRow, 0, upperLimit + maxDiff + 1)
      ),
      fill = '#f0f0f0',
      stat = 'identity'
    )
  
  for (y in seq(1, maxDiff)) {
    plot = plot + geom_hline(yintercept = y, linetype = 'dashed', colour = '#808080')
    plot = plot + geom_hline(yintercept = -y, linetype = 'dashed', colour = '#808080')
  }
  
  plot = plot +
    geom_hline(yintercept = 0, colour = '808080') +
    stat_summary(
      mapping = aes(y = `Vela => HIVDB`),
      geom="boxplot",
      width=0.4,
      fun.data=function(x) {
        r = quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
        names(r) = c("ymin", "lower", "middle", "upper", "ymax")
        r
      }
    ) +
    geom_text(
      mapping = aes(label = num_unusual),
      y = offsetNumUnusual + 1.1,
      hjust = 0.5,
      size = 3
    ) +
    geom_text(
      mapping = aes(label = num_apobec),
      y = offsetNumApobec + 1.1,
      hjust = 0.5,
      size = 3
    ) +
    geom_text(
      mapping = aes(label = num_stop),
      y = offsetNumStop + 1.1,
      hjust = 0.5,
      size = 3
    ) +
    geom_text(
      mapping = aes(label = ifelse(is.na(VelaMuts), "", VelaMuts)),
      y = offsetVelaOnly,
      hjust = 0,
      size = 3,
      colour = "#d93949"
    ) +
    geom_text(
      mapping = aes(label = ifelse(is.na(HIVDBMuts), "", HIVDBMuts)),
      y = offsetHIVDBOnly,
      hjust = 0,
      size = 3,
      colour = "#25a545"
    ) +
    scale_x_discrete(
      name = "Sample",
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "Vela => HIVDB",
      expand = c(0, 0),
      limits = c(-maxDiff - 0.5, upperLimit + 0.5),
      breaks=seq(-maxDiff, maxDiff),
      sec.axis = sec_axis(
        trans = ~., name = "Vela => HIVDB",
        breaks = c(seq(-maxDiff, maxDiff),
                   offsetNumUnusual + 1.1,
                   offsetNumApobec + 1.1,
                   offsetNumStop + 1.1,
                   offsetVelaOnly + 0.9,
                   offsetHIVDBOnly + 1.1),
        labels = c(seq(-maxDiff, maxDiff),
                   "# Unusuals",
                   "# APOBECs",
                   "# Stops",
                   "Vela DRMs",
                   "HIVDB DRMs")
      )
    ) +
    coord_flip(clip = "off") +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  plot
}

plot.correlationAgainstLevelDiffs = function(data, x, hide_y_axis_title = FALSE) {
  x = eval(substitute(x), data, parent.frame())
  plot = ggplot(data, aes(x = x, y = `Vela => HIVDB`)) +
    geom_count() + scale_size_area() +
    scale_y_continuous(
      name = 'Level diff',
      limits = c(-4, 4),
      breaks = seq(-4, 4),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      breaks = seq(
        floor(min(x) / 5) * 5,
        max(x),
        max(1, floor((max(x) - min(x)) / 8))
      )
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.background =element_blank(),
      panel.border = element_rect(fill = 'transparent', colour = "black")
    )
  if (hide_y_axis_title) {
    plot = plot + theme(axis.title.y = element_blank())
  }
  plot
}