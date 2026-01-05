library(ggplot2)
library(cowplot)

Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "Data")


RunNames <- list("NoSpaceBinaryInterac", "BinaryInterac", "ContinuousInterac", "RareSpecies")

plotlist <- lapply(RunNames, FUN = function(RunName) {
    ## Data
    load(file.path(Dir.Data, paste0(RunName, "_DEMO_Environment.RData")))
    load(file.path(Dir.Data, paste0(RunName, "_DEMO.RData")))

    ## Environment
    Env_long <- as.data.frame(
        as.table(
            Env_mat[
                seq(from = 1, to = nrow(Env_mat), by = 9),
                seq(from = 1, to = ncol(Env_mat), by = 9)
            ]
        )
    )
    colnames(Env_long) <- c("Y", "X", "VALUES")
    Env_long <- apply(Env_long, 2, as.numeric)

    Env_gg <- ggplot(Env_long, aes(x = X, y = Y, fill = VALUES)) +
        geom_tile() +
        coord_fixed() +
        scale_fill_viridis_c(
            option = "C",
            name = "Environmental Value /\n Optimal Phenotype",
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        theme_bw() +
        theme(
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            legend.position = "top",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c")
        )
    # Env_gg

    ## Network
    Effect_Mat <- igraph::as_adjacency_matrix(Network_igraph, attr = "weight")
    rownames(Effect_Mat) <- colnames(Effect_Mat) <- names(CarryingK_vec)
    edg_df <- as.data.frame(as.table(as.matrix(Effect_Mat)))
    colnames(edg_df) <- c("Partner 1", "Partner 2", "Strength")
    NetMat_gg <- ggplot(edg_df, aes(x = `Partner 1`, y = `Partner 2`, fill = Strength)) +
        geom_tile(color = "black", lwd = 0.5, linetype = 1) +
        coord_fixed() +
        scale_fill_gradient2(
            low = "#5ab4ac",
            high = "#d8b365",
            name = "Association Strength",
            guide = guide_colourbar(title.vjust = 0.75)
        ) +
        # scale_fill_viridis_c(
        #     option = "C",
        #     name = "Environmental Value /\n Optimal Phenotype",
        #     guide = guide_colourbar(title.vjust = 0.75)
        # ) +
        theme_bw() +
        theme(
            plot.margin = unit(c(0, 0, 0, 0), "cm"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.width = unit(2, "cm"),
            legend.key.height = unit(1, "cm"),
            panel.background = element_rect(fill = "#2c2c2c", color = "#2c2c2c"),
            axis.text.x = element_text(angle = -20, hjust = 0)
        )
    # NetMat_gg

    ## Final spatial arrangement
    last_df <- SimResult[[length(SimResult)]]
    last_df$Species <- factor(as.numeric(gsub(last_df$Species, pattern = "Sp_", replacement = "")))
    Traits_df <- aggregate(Trait ~ Species, data = last_df, FUN = mean)
    Traits_df$Species <- factor(as.numeric(gsub(Traits_df$Species,
        pattern = "Sp_", replacement = ""
    )))
    Final_gg <- ggplot(
        last_df,
        aes(x = X, y = Y, col = Species, shape = Species)
    ) +
        geom_point() +
        scale_shape_manual(values = 1:nlevels(last_df$Species)) +
        scale_color_viridis_d() +
        geom_vline(
            data = Traits_df,
            aes(
                xintercept = Trait,
                col = Species
            )
        ) +
        xlim(
            0,
            10
        ) +
        theme_bw() +
        theme(legend.position = "top")
    # Final_gg

    ## Abundance over time
    Abund_time <- lapply(names(SimResult),
        FUN = function(t) {
            ID_iter <- SimResult[[t]]
            cbind(data.frame(table(ID_iter$Species)), t)
        }
    )
    Abund_time <- do.call(rbind, Abund_time)
    Abund_time$t <- as.numeric(Abund_time$t)
    colnames(Abund_time)[1:2] <- c("Species", "Abundance")
    Abund_time$Species <- factor(as.numeric(gsub(Abund_time$Species,
        pattern = "Sp_", replacement = ""
    )))

    AbundTime_gg <- ggplot(Abund_time, aes(x = t, y = Abundance, col = Species)) +
        geom_line(linewidth = 1.5) +
        scale_color_viridis_d() +
        theme_bw() +
        theme(legend.position = "bottom")
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
    file = "DemoFig.png",
    width = 28, height = 34
)
