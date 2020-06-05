# Two-way dotplot function
require(Seurat)
require(cowplot)

DotPlot2 <- function(
  object,
  features,
  group1,
  group2,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  scale = TRUE,
  scale.min = NA,
  scale.max = NA
) {
  data.features <- FetchData(object = object, vars = features)
  data.features$grp1 <- object[[group1, drop = TRUE]]
  data.features$grp2 <- object[[group2, drop = TRUE]]
  if (!is.factor(x = data.features$grp1)) {
    data.features$grp1 <- factor(x = data.features$grp1)
  }
  if (!is.factor(x = data.features$grp2)) {
    data.features$grp2 <- factor(x = data.features$grp2)
  }
  grp1.levels <- levels(x = data.features$grp1)
  grp2.levels <- levels(x = data.features$grp2)
  data.features$grp1 <- as.vector(x = data.features$grp1)
  data.features$grp2 <- as.vector(x = data.features$grp2)
  
  data.plot <- lapply(grp1.levels,
    function(grp1) {
      lapply(grp2.levels, function(grp2) {
        data.use <- data.features[data.features$grp1 == grp1 & data.features$grp2 == grp2, 1:(ncol(x = data.features) - 2), drop = FALSE]
        avg.exp <- apply(
          X = data.use,
          MARGIN = 2,
          FUN = function(x) {
            return(mean(x = expm1(x = x)))
          })
        pct.exp <- apply(data.use, 2, function(x) length(x[x > 0]) / length(x = x))
        return(list(avg.exp = avg.exp, pct.exp = pct.exp))
      })
    }
  )
  data.plot = setNames(lapply(data.plot, setNames, grp2.levels), grp1.levels) %>% unlist(recursive = F)
  
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled.")
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log(x = data.use)
      }
      return(data.use)
    }
  )
  
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = rev(x = features)
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  color.by <- 'avg.exp.scaled'
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  data.plot$grp1 = sapply(data.plot$id, function(x) strsplit(x, split = "\\.")[[1]][1]) %>% factor(levels = grp1.levels)
  data.plot$grp2 = sapply(data.plot$id, function(x) strsplit(x, split = "\\.")[[1]][2]) %>% factor(levels = grp2.levels)
  plots <- lapply(features, 
                  function(x) {
                    data.use = data.plot[data.plot$features.plot == x,]
                    ggplot(data = data.use, mapping = aes_string(x = 'grp1', y = 'grp2')) +
                    geom_point(mapping = aes_string(size = 'pct.exp', color = color.by)) +
                    scale_color_gradient(low = cols[1], high = cols[2], limits = c(col.min, col.max)) +
                    scale_radius(range = c(0, dot.scale), limits = c(0, 100)) +
                    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
                    guides(size = guide_legend(title = 'Percent Expressed'),
                           color = guide_colorbar(title = 'Average Expression')) +
                    labs(title = x, x = group1,y = group2) +
                    theme_cowplot()
                  })
  return(plots)
}