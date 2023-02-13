suppressMessages(library(dplyr))
suppressMessages(library(formattable))

fullTableOfStats<-read.table(file = snakemake@input[[1]], sep = "", header = TRUE, row.names=NULL, check.names = FALSE)

fullTableOfStats$Gaps_per_Gb <- ((fullTableOfStats$Gaps * fullTableOfStats$Est_Size_Mb)/1000)
fullTableOfStats$Gaps_per_Gb <- as.integer(fullTableOfStats$Gaps_per_Gb)



selectionOfStats_internalComparisonHeatmap <- fullTableOfStats %>%
  select(c('ASM_ID','ASM_LEVEL', 'Bases_Mb', 'Het_%','GC_%', 'Gaps_per_Gb', 'Scaff', 'Cont',
           'Longest_Scaff_Mb','Scaff_NG50_Mb', 'Scaff_NG95_Mb',
           'Longest_Cont_Mb', 'Cont_NG50_Mb','Cont_NG95_Mb',
           'QV', 'Completeness',
           'Comp_BUSCOs_%','Comp_Single_BUSCOs_%'))




customBlue_max = "#2e96ff"

customBlue_min = "#dcecfc"


customGray_max = "#8c8c8c"

customGray_min = "#e6e6e6"

sink(file = snakemake@output[[1]])

format_table(selectionOfStats_internalComparisonHeatmap,
             align =c("l","c","c","c","c", "c", "c", "c", "c"),
             list(Bases_Mb = color_tile(customGray_min, customGray_max),
                  `Het_%` = color_tile(customGray_min, customGray_max),
                  `GC_%` = color_tile(customGray_min, customGray_max),
                  Gaps_per_Gb = color_tile(customBlue_max, customBlue_min),
                  Scaff = color_tile(customBlue_max,customBlue_min),
                  Cont = color_tile(customBlue_max,customBlue_min),
                  Longest_Scaff_Mb = color_tile(customBlue_min, customBlue_max),
                  Scaff_NG50_Mb = color_tile(customBlue_min, customBlue_max),
                  Scaff_NG95_Mb = color_tile(customBlue_min, customBlue_max),
                  Longest_Cont_Mb = color_tile(customBlue_min, customBlue_max),
                  Cont_NG50_Mb = color_tile(customBlue_min, customBlue_max),
                  Cont_NG95_Mb = color_tile(customBlue_min, customBlue_max),
                  QV = color_tile(customBlue_min, customBlue_max),
                  Completeness = color_tile(customBlue_min, customBlue_max),
                  `Comp_BUSCOs_%` = color_tile(customBlue_min, customBlue_max),
                  `Comp_Single_BUSCOs_%` = color_tile(customBlue_min, customBlue_max)))

cat('<br>')
cat("\n")
cat('<br>')

legendTable<-read.table(file = snakemake@params[[1]], sep = "", header = TRUE, row.names=NULL, check.names = FALSE)

format_table(legendTable,
            align =c("c","c","c", "c", "c", "c", "c","c","c","c", "c", "c", "c", "c","c", "c"),
            list(Bases_Mb = formatter("span",
                                         style = ~style(display = "block",
                                                        padding = "0 4px",
                                                        `border-radius` = "3px",
                                                        `background-color` = ifelse(Bases_Mb == 'Max',customGray_max, customGray_min))),
                 `Het_%` = formatter("span",
                                      style = ~style(display = "block",
                                                     padding = "0 4px",
                                                     `border-radius` = "3px",
                                                     `background-color` = ifelse(`Het_%` == 'Max',customGray_max, customGray_min))),
                 `GC_%` = formatter("span",
                                     style = ~style(display = "block",
                                                    padding = "0 4px",
                                                    `border-radius` = "3px",
                                                    `background-color` = ifelse(`GC_%` == 'Max',customGray_max, customGray_min))),
                 Gaps_per_Gb = formatter("span",
                                    style = ~style(display = "block",
                                                   padding = "0 4px",
                                                   `border-radius` = "3px",
                                                   `background-color` = ifelse(Gaps_per_Gb == 'Min',customBlue_max, customBlue_min))),
                 Scaff = formatter("span",
                                         style = ~style(display = "block",
                                                        padding = "0 4px",
                                                        `border-radius` = "3px",
                                                        `background-color` = ifelse(Scaff == 'Min',customBlue_max, customBlue_min))),
                 Cont = formatter("span",
                                         style = ~style(display = "block",
                                                        padding = "0 4px",
                                                        `border-radius` = "3px",
                                                        `background-color` = ifelse(Cont == 'Min',customBlue_max, customBlue_min))),
                 Longest_Scaff_Mb = formatter("span",
                                    style = ~style(display = "block",
                                                 padding = "0 4px",
                                                 `border-radius` = "3px",
                                                 `background-color` = ifelse(Longest_Scaff_Mb == 'Max',customBlue_max, customBlue_min))),
                 Scaff_NG50_Mb = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(Scaff_NG50_Mb == 'Max',customBlue_max, customBlue_min))),
                 Scaff_NG95_Mb = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(Scaff_NG95_Mb == 'Max',customBlue_max, customBlue_min))),
                 Longest_Cont_Mb = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(Longest_Cont_Mb == 'Max',customBlue_max, customBlue_min))),
                 Cont_NG50_Mb = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(Cont_NG50_Mb == 'Max',customBlue_max, customBlue_min))),
                 Cont_NG95_Mb = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(Cont_NG95_Mb == 'Max',customBlue_max, customBlue_min))),
                 QV = formatter("span",
                                              style = ~style(display = "block",
                                                             padding = "0 4px",
                                                             `border-radius` = "3px",
                                                             `background-color` = ifelse(QV == 'Max',customBlue_max, customBlue_min))),
                 Completeness = formatter("span",
                                style = ~style(display = "block",
                                               padding = "0 4px",
                                               `border-radius` = "3px",
                                               `background-color` = ifelse(Completeness == 'Max',customBlue_max, customBlue_min))),
                 `Comp_BUSCOs_%` = formatter("span",
                                style = ~style(display = "block",
                                               padding = "0 4px",
                                               `border-radius` = "3px",
                                               `background-color` = ifelse(`Comp_BUSCOs_%` == 'Max',customBlue_max, customBlue_min))),
                 `Comp_Single_BUSCOs_%` = formatter("span",
                                style = ~style(display = "block",
                                               padding = "0 4px",
                                               `border-radius` = "3px",
                                               `background-color` = ifelse(`Comp_Single_BUSCOs_%` == 'Max',customBlue_max, customBlue_min)))))

sink(file = NULL)
