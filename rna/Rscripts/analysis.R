vicar_l2c_results <- read_csv(here::here("rna/output/crvi-ctrl-DE_vicar-table_2023-12-04.csv"))


# DE state -- upregulated, downregulated, or not DE --
vicar_l2c_results <- vicar_l2c_results %>% 
    mutate(sig = case_when(
        betahat > 0 & lfdr < 0.05 ~ "up",
        betahat < 0 & lfdr < 0.05 ~ "down",
        TRUE ~ "not"
    ))

# Get overlaps between significant genes between tissues
# (i.e. genes that are significant in both tissues)

sorted_vicar_l2c_results <- vicar_l2c_results %>% 
    arrange(lfdr)

body_sig <- vicar_l2c_results %>% 
    filter(tissue == "body", sig != "not")

head_sig <- vicar_l2c_results %>%
    filter(tissue == "head", sig != "not") %>%
    arrange(lfdr)

overlap_genes <- intersect(body_sig$Geneid, head_sig$Geneid)


# get only genes that are significant in both tissues
vicar_l2c_overlap <- vicar_l2c_results %>% 
    filter(Geneid %in% overlap_genes)


# Plot lfdr against effect size --
volcano <- ggplot(vicar_l2c_results, aes(x = betahat, y = -log10(lfdr), color = sig, label =Geneid  )) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_bw() +
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), # to set the colours of our variable<br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
    labs(x = "Estimated Effect Size", y = "-log10(Local FDR)") +
    facet_wrap(~tissue, ncol = 2, scales = "free_y") +
    geom_text_repel(data = subset(vicar_l2c_results, lfdr < 0.05), 
                    aes(label = Geneid),
                    max.overlaps = 5,
                    size = 3, 
                    nudge_x = 0.5, 
                    nudge_y = 0.5,
                    segment.alpha = 0.2,
                    segment.color = "grey50",
                    seed = 123,
                    show.legend = FALSE)

save_plot(here::here(append_date("rna/output/crvi-ctrl-DE_vicar-lfdr_vs_logFC.png")),volcano, base_width = 6, base_height = 5, ncol = 2,dpi=300)

x = vicar_l2c_results |>
    dlply("tissue") |> 
    map(\(x) slice_min(x, lfdr, n = 100)) |> 
    map(\(x) mutate(x, n = 1:nrow(x))) |>
    list_rbind() |>
    select(n, Geneid, tissue, lfdr)

p = ggplot(x, aes(n, y = lfdr)) + 
    geom_point() + 
    facet_wrap(~tissue, scales = "free_y") +
    theme_bw()

save_plot("tmp/test.png",p , base_width = 7, base_height = 7, ncol = 2)