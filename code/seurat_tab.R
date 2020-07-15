# Create table from seurat object
# Date: 05/25/2020
# Author: Mike Cuoco

library(tidyverse)
library(tidyselect)
library(reshape2)
library(gt)

seurat_tab <- function(seurat_object, col_var, row_var, group_var=NULL, title){
  # Sort input columns, only works with numeric vector currently, needs to inherit levels of parent
  cols = as.character(sort(unique(seurat_object[[]][[col_var]])))
  
  if (is.null(group_var)){
    form = formula(paste0("total + ",row_var," ~ ",col_var))
    tab = seurat_object[[]] %>%
      group_by(!!sym(col_var), !!sym(row_var)) %>% 
      tally() %>%
      ungroup() %>%
      ddply(.(eval(sym(row_var))), mutate, total = sum(n)) %>%
      dcast(form, value.var = "n", sum, fill=0) %>%
      arrange(eval(sym(row_var))) %>%
      gt(rowname_col = row_var) %>%
      tab_header(title = title, subtitle = paste("dataset:",seurat_object@project.name)) %>%
      tab_spanner(label = col_var, columns = tidyselect::any_of(cols)) %>%
      cols_move_to_end("total") %>%
      summary_rows(fns = list(total = "sum"), decimals = 0) %>%
      tab_options(
        summary_row.background.color = "#ACEACE",
        grand_summary_row.background.color = "#ACEACE",
        heading.title.font.size = 25,
        heading.title.font.weight = "bold",
        column_labels.background.color = "#B4DCF7",
        column_labels.font.size = 20,
        column_labels.font.weight = "bold"
      )
  }
  else {
    form = formula(paste0("total + ",group_var," + ",row_var," ~ ",col_var))
    tab = seurat_object[[]] %>%
            group_by(!!sym(col_var), !!sym(row_var), !!sym(group_var)) %>% 
            tally() %>%
            ungroup() %>%
            ddply(.(eval(sym(row_var)), eval(sym(group_var))), mutate, total = sum(n)) %>%
            dcast(form, value.var = "n", sum, fill=0) %>%
            group_by(!!sym(group_var)) %>%
            arrange(eval(sym(row_var))) %>%
            gt(rowname_col = row_var, groupname_col = group_var) %>%
            tab_header(title = title, subtitle =  paste("dataset:",seurat_object@project.name)) %>%
            tab_spanner(label = col_var, columns = tidyselect::any_of(cols)) %>%
            cols_move_to_end("total") %>%
            summary_rows(groups = T, fns = list(total = "sum"), decimals = 0) %>%
            grand_summary_rows(fns = list(`grand total` = "sum"), decimals = 0) %>%
            tab_options(
              summary_row.background.color = "#ACEACE",
              grand_summary_row.background.color = "#ACEACE",
              row_group.background.color = "#FFEFDB",
              heading.title.font.size = 25,
              heading.title.font.weight = "bold",
              column_labels.background.color = "#B4DCF7",
              column_labels.font.size = 20,
              column_labels.font.weight = "bold"
            )
  }
  return(tab)
}
