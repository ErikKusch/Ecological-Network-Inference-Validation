plotlist <- lapply(RunNames, FUN = function(RunName) {
    # RunName <- RunNames[1]
    # print(RunName)
    ## Data
    load(file.path(Dir.Data, paste0(RunName, "_DEMO_Environment.RData")))
    load(file.path(Dir.Data, paste0(RunName, "_DEMO.RData")))

    ## Environment
    Env_gg <- Plot_Environment(Env_mat)
    # Env_gg

    ## Network
    NetMat_gg <- Plot_NetMat(Network_igraph)
    # NetMat_gg

    ## Final spatial arrangement
    last_df <- SimResult[[length(SimResult)]]
    Final_gg <- Plot_IndivsInSpace(last_df)
    # Final_gg

    ## Abundance over time
    AbundTime_gg <- Plot_AbundTime_gg(SimResult)
    # AbundTime_gg

    ## return
    list(
        Env = Env_gg,
        Net = NetMat_gg,
        Placement = Final_gg,
        Abund_time = AbundTime_gg
    )
})
names(plotlist) <- unlist(RunNames)

## fuse plots
label_row <- function(text) {
    ggplot() +
        geom_rect(
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = "white",
            color = NA
        ) +
        annotate(
            "text",
            x = 0, y = 0,
            label = text,
            hjust = 0,
            fontface = "bold",
            size = 5
        ) +
        coord_cartesian(xlim = c(0, 1), ylim = c(-1, 1), clip = "off") +
        theme_void() +
        theme(
            plot.margin = margin(0, 0, 0, 0)
        )
}

label_A <- label_row("(A) Environment")
label_B <- label_row("(B) Network Matrix")
label_C <- label_row("(C) Final Placement of Individuals in Environment")
label_D <- label_row("(D) Abundance of Species throughout the Simulation")

main_ggs <- plot_grid(
    label_A,
    plot_grid(plotlist = lapply(plotlist, "[[", "Env"), nrow = 1),
    label_B,
    plot_grid(plotlist = lapply(plotlist, "[[", "Net"), nrow = 1),
    label_C,
    plot_grid(plotlist = lapply(plotlist, "[[", "Placement"), nrow = 1),
    label_D,
    plot_grid(plotlist = lapply(plotlist, "[[", "Abund_time"), nrow = 1),
    ncol = 1,
    rel_heights = c(
        0.08, 1, # A label + plot
        0.08, 1, # B label + plot
        0.08, 1, # C label + plot
        0.08, 1 # D label + plot
    )
)
main_ggs
ggsave(
    main_ggs,
    file = file.path(Dir.Exports, "Figure2_SimulationRunsDemo.png"),
    width = 28, height = 34
)
