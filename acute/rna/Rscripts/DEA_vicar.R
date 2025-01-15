table <- 
     list_rbind(list(head = mutate(rnaseq_data$head$mwash$l2c$result,
                                   Geneid = rownames(rnaseq_data$head$mwash$l2c$result)),
                    body = mutate(rnaseq_data$body$mwash$l2c$result,
                                   Geneid = rownames(rnaseq_data$body$mwash$l2c$result))), 
               names_to = "tissue") %>% 
     as_tibble %>% 
     relocate(Geneid)

write_csv(table, affix_date("output/crvi-ctrl-DE_vicar-table.csv"))
