suppressMessages(library(dplyr))
suppressMessages(library(formattable))
#suppressMessages(library(data.table))


fullTableOfStats<-read.table(file = snakemake@input[[1]], sep = "", header = TRUE, row.names=NULL, check.names = FALSE)



fullTableOfStats$Gaps_per_Gb <- ((fullTableOfStats$Gaps * fullTableOfStats$Est_Size_Mb)/1000)
fullTableOfStats$Gaps_per_Gb <- as.integer(fullTableOfStats$Gaps_per_Gb)


customPurple='#c699e8'
customGreen='#8fc773'

selectionOfStats_colouredHeatmap <- fullTableOfStats %>%
  select(c('ASM_ID','ASM_LEVEL',
           'Gaps_per_Gb', 'Scaff_NG50_Mb',
           'Cont_NG50_Mb','QV',
           'Completeness', 'Comp_Single_BUSCOs_%'))

sink(file = snakemake@output[[1]])
format_table(selectionOfStats_colouredHeatmap,
             align =c("l","c","c","c","c", "c", "c", "c", "c"),
             list(Gaps_per_Gb = formatter("span",
                                          			style = ~style(display = "block",
                                                         	padding = "0 4px",
                                                         	`border-radius` = "3px",
                                                         	`background-color` = ifelse(ASM_LEVEL == 'cont', "#666666",
                                                                                     	ifelse(between(Gaps_per_Gb,1001, 10000), "#FFCC99",
                                                                                            ifelse(between(Gaps_per_Gb,201 , 1000), customGreen,
                                                                                                   ifelse(Gaps_per_Gb <= 200, customPurple, "#FF9999")))))),
                  Scaff_NG50_Mb = formatter("span",
                                            			style = ~style(display = "block",
                                                           	padding = "0 4px",
                                                           	`border-radius` = "3px",
                                                           	`background-color` = ifelse(between(Scaff_NG50_Mb,0.1, 9.99999), "#FFCC99",
                                                                                       ifelse(between(Scaff_NG50_Mb,10, 99.99999), customGreen,
                                                                                              ifelse(Scaff_NG50_Mb >= 100, customPurple, "#FF9999"))))),
                  Cont_NG50_Mb = formatter("span",
                                           			style = ~style(display = "block",
                                                          	padding = "0 4px",
                                                          	`border-radius` = "3px",
                                                          	`background-color` = ifelse(between(Cont_NG50_Mb,0.01, 0.99999), "#FFCC99",
                                                                                      ifelse(between(Cont_NG50_Mb,1, 9.999999), customGreen,
                                                                                             ifelse(Cont_NG50_Mb >= 10, customPurple, "#FF9999"))))),
                  Completeness = formatter("span",
                                           			style = ~style(display = "block",
                                                          	padding = "0 4px",
                                                          	`border-radius` = "3px",
                                                          	`background-color` = ifelse(between(Completeness,80, 89.9999), "#FFCC99",
                                                                                      ifelse(between(Completeness,90, 94.99999), customGreen,
                                                                                             ifelse(Completeness >= 95, customPurple, "#FF9999"))))),
                  `Comp_Single_BUSCOs_%` = formatter("span",
                                                 		style = ~style(display = "block",
                                                                padding = "0 4px",
                                                                `border-radius` = "3px",
                                                                `background-color` = ifelse(between(`Comp_Single_BUSCOs_%`,80, 89.9999), "#FFCC99",
                                                                                            ifelse(between(`Comp_Single_BUSCOs_%`,90, 95), customGreen,
                                                                                                   ifelse(`Comp_Single_BUSCOs_%` >= 95, customPurple, "#FF9999"))))),
                  QV = formatter("span",
                                 				style = ~style(display = "block",
                                                		padding = "0 4px",
                                                		`border-radius` = "3px",
                                                		`background-color` = ifelse(between(QV,35, 39.9999999999), "#FFCC99",
                                                                            		ifelse(between(QV,40, 49.999999999), customGreen,
                                                                                   		ifelse(QV >= 50, customPurple, "#FF9999")))))))

cat('<br>')
cat("\n")
cat('<br>')

legendTable<-read.table(file = snakemake@params[[1]], sep = "", header = TRUE, row.names=NULL, check.names = FALSE)

format_table(legendTable,
            align =c("c","c","c", "c", "c", "c", "c"),
            list(Gaps_per_Gb = formatter("span",
                                         style = ~style(display = "block",
                                                        padding = "0 4px",
                                                        `border-radius` = "3px",
							"font.size" = "12px",
 							"width" = '150px',
                                                        `background-color` = ifelse(Gaps_per_Gb == '> 10000',"#FF9999",
                                                                                    ifelse(Gaps_per_Gb == '1000 - 10000',"#FFCC99",
                                                                                           ifelse(Gaps_per_Gb == '200 - 1000',customGreen, customPurple))))),
                 Scaff_NG50_Mb = formatter("span",
                                         style = ~style(display = "block",
                                                        padding = "0 4px",
                                                        `border-radius` = "3px",
							"font.size" = "12px",
							"width" = '150px',
                                                        `background-color` = ifelse(Scaff_NG50_Mb == '< 0.1Mbp',"#FF9999",
                                                                                    ifelse(Scaff_NG50_Mb == '0.1Mbp - 10Mbp',"#FFCC99",
                                                                                           ifelse(Scaff_NG50_Mb == '10Mbp - 100Mbp',customGreen, customPurple))))),
                 Cont_NG50_Mb = formatter("span",
                                           style = ~style(display = "block",
                                                          padding = "0 4px",
                                                          `border-radius` = "3px",
							  "font.size" = "12px",
							  "width" = '150px',
                                                          `background-color` = ifelse(Cont_NG50_Mb == '< 0.01Mbp',"#FF9999",
                                                                                      ifelse(Cont_NG50_Mb == '0.01Mbp - 1Mbp',"#FFCC99",
                                                                                             ifelse(Cont_NG50_Mb == '1Mbp - 10Mbp',customGreen, customPurple))))),
                 QV = formatter("span",
                                          style = ~style(display = "block",
                                                         padding = "0 4px",
                                                         `border-radius` = "3px",
							 "font.size" = "12px",
							 "width" = '100px',
                                                         `background-color` = ifelse(QV == '< 35',"#FF9999",
                                                                                     ifelse(QV == '35 - 40',"#FFCC99",
                                                                                            ifelse(QV == '40 - 50',customGreen, customPurple))))),
                 Completeness = formatter("span",
                                style = ~style(display = "block",
                                               padding = "0 4px",
                                               `border-radius` = "3px",
					       "font.size" = "12px",
					       "width" = '150px',
                                               `background-color` = ifelse(Completeness == '< 80%',"#FF9999",
                                                                           ifelse(Completeness == '80% - 90%',"#FFCC99",
                                                                                  ifelse(Completeness == '90% - 95%',customGreen, customPurple))))),
                 `Comp_Single_BUSCOs_%` = formatter("span",
                                          style = ~style(display = "block",
                                                         padding = "0 4px",
                                                         `border-radius` = "3px",
							 "font.size" = "12px",
							 "width" = '220px',
                                                         `background-color` = ifelse(`Comp_Single_BUSCOs_%` == '< 80%',"#FF9999",
                                                                                     ifelse(`Comp_Single_BUSCOs_%` == '80% - 90%',"#FFCC99",
                                                                                            ifelse(`Comp_Single_BUSCOs_%` == '90% - 95%',customGreen, customPurple)))))))
sink(file = NULL)
