# NO NETWORK SPECIES IN ENVIRONMENT =========================================
Network_igraph <- Sim.Network(
    n_spec = 4,
    NetworkType = "Association",
    Sparcity = 1,
    MaxStrength = 0, # what absolute value is the maximum link strength
    seed = 42 # seed for randomness
)

CarryingK_vec <- c(Sp_1 = 300, Sp_2 = 300, Sp_3 = 300, Sp_4 = 300)

Niches_vec <- c(Sp_1 = 2, Sp_2 = 4, Sp_3 = 6, Sp_4 = 8)

Initialise_df <- Sim.Initialise(
    n_spec = 4,
    n_individuals = 1e3,
    n_mode = "each",
    Env_range = c(0, 10),
    Trait_means = Niches_vec,
    Trait_sd = 0,
    seed = 42
)
# Plot_IndivsInSpace(Initialise_df)

Env_mat <- Sim.Space(
    x_range = c(0, 10),
    y_range = c(0, 10),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) x,
    y_gradient = function(y) y
) / 2
# Plot_Environment(Env_mat)

FNAME <- file.path(Dir.Concept, "BoxDEMO_SpaceOnly_SimResult.RData")

if (file.exists(FNAME)) {
    load(FNAME)
} else {
    SimResult <- Sim.Compute(
        # Demographic parameters
        d0 = d0,
        b0 = b0,
        k_vec = CarryingK_vec,
        ID_df = Initialise_df,

        # Spatial parameters
        env.xy = Env_mat,
        env.sd = 0.5,
        mig.sd = 0.5,
        mig.top = 0.35,

        # Interaction parameters
        interac.maxdis = 1,
        interac.igraph = Network_igraph,
        interac.scale = 0,

        # Simulation parameters
        Sim.t.max = t_max * 2,
        Sim.t.inter = t_inter,
        seed = 0,
        verbose = verbose, # whether to print progress in time as current time
        RunName = "BoxDEMO_SpaceOnly"
    )
    save(SimResult, file = FNAME)
}
last_df <- SimResult[[length(SimResult)]]

p_main <- Plot_Environment(Env_mat) +
    geom_point(
        data = last_df,
        aes(x = X, y = Y, shape = Species),
        fill = "white"
    ) +
    scale_shape_manual(values = 1:length(unique(last_df$Species)))

p_shape <- ggplot(
    last_df,
    aes(x = X, y = Y, shape = Species)
) +
    geom_point() +
    scale_shape_manual(values = 1:length(unique(last_df$Species))) +
    theme_void() +
    theme(legend.position = "right")

shape_legend <- get_legend(p_shape)

p_main_clean <- p_main +
    guides(shape = "none")

Enviro_SpecinSpace_gg <- plot_grid(
    p_main_clean,
    shape_legend,
    ncol = 2,
    rel_widths = c(1, 0.1)
)
ggsave(
    filename = file.path(Dir.Exports, "Figure_BoxDEMO_SpaceOnly_SpeciesInSpace.jpg"),
    plot = Enviro_SpecinSpace_gg,
    width = 18, height = 18, units = "cm"
)
# print(Enviro_SpecinSpace_gg)

# 1-1 ASSOCIATIONS WITH ENVIRONMENT ========================================
AssocPanels <- list("Competition" = -1, "Mutualism" = 1, "No Association" = 0)

AssocPanels_gg <- lapply(1:length(AssocPanels), function(x) {
    x <- AssocPanels[[x]]
    if (x == 1) {
        y <- c("Normal", "OvercomingEnvironment")
    } else {
        y <- "Normal"
    }
    # print(x)

    plots <- lapply(y, function(i) {
        FNAME <- file.path(Dir.Concept, paste0("BoxDEMO_AssocsPanels_SimResult_Assoc_", x, gsub(i, pattern = "Normal", replacement = ""), ".RData"))
        # print(FNAME)

        Network_igraph <- Sim.Network(
            n_spec = 2,
            NetworkType = "Association",
            Sparcity = 0,
            MaxStrength = 1, # what absolute value is the maximum link strength
            seed = 42 # seed for randomness
        )
        E(Network_igraph)$weight <- x

        CarryingK_vec <- rep(200, length(V(Network_igraph)))
        names(CarryingK_vec) <- paste0("Sp_", V(Network_igraph))

        if (i == "OvercomingEnvironment") {
            Niches_vec <- c(Sp_1 = 1.5, Sp_2 = 3.5)
        } else {
            Niches_vec <- c(Sp_1 = 0, Sp_2 = 0)
        }

        Initialise_df <- Sim.Initialise(
            n_spec = 2,
            n_individuals = 2e3,
            n_mode = "each",
            Env_range = c(0, 10),
            Trait_means = Niches_vec,
            Trait_sd = 0,
            seed = 42
        )
        # Plot_IndivsInSpace(Initialise_df)

        if (i == "OvercomingEnvironment") {
            Env_mat <- Sim.Space(
                x_range = c(0, 10),
                y_range = c(0, 10),
                ncol = 1e3, nrow = 1e3,
                x_gradient = function(x) x,
                y_gradient = function(y) 0
            ) / 2
        } else {
            Env_mat <- Sim.Space(
                x_range = c(0, 10),
                y_range = c(0, 10),
                ncol = 1e3, nrow = 1e3,
                x_gradient = function(x, alpha = 1, L = 10) {
                    exp(-alpha * x) + exp(-alpha * (L - x))
                },
                y_gradient = function(y, alpha = 1, L = 10) {
                    exp(-alpha * y) + exp(-alpha * (L - y))
                }
            )
        }
        # Plot_Environment(Env_mat)

        if (file.exists(FNAME)) {
            load(FNAME)
        } else {
            SimResult <- Sim.Compute(
                # Demographic parameters
                d0 = d0,
                b0 = b0,
                k_vec = CarryingK_vec,
                ID_df = Initialise_df,

                # Spatial parameters
                env.xy = Env_mat,
                env.sd = 0.5,
                mig.sd = 0.25,
                mig.top = 0.15,

                # Interaction parameters
                interac.maxdis = 0.5,
                interac.igraph = Network_igraph,
                interac.scale = 1,

                # Simulation parameters
                Sim.t.max = t_max,
                Sim.t.inter = t_inter,
                seed = 42,
                verbose = verbose, # whether to print progress in time as current time
                RunName = "BoxDEMO_AssocsPanels"
            )
            save(SimResult, file = FNAME)
        }

        last_df <- SimResult[[length(SimResult)]]

        p_main <- Plot_Environment(Env_mat) +
            geom_point(
                data = last_df,
                aes(x = X, y = Y, shape = Species),
                fill = "white", col = "white"
            ) +
            scale_shape_manual(values = c(3, 5)) +
            theme(legend.position = "top")
        p_main
    })
    plots
})
names(AssocPanels_gg) <- names(AssocPanels)

legend <- get_legend(AssocPanels_gg[[1]][[1]])

AssocPanels_gg <- unlist(AssocPanels_gg[c(3, 1, 2)])

AssocPanels_gg2 <- plot_grid(
    legend,
    plot_grid(plotlist = lapply(AssocPanels_gg, FUN = function(x) {
        x + theme(legend.position = "none") + labs(x = "", y = "")
    }), ncol = 2),
    ncol = 1,
    rel_heights = c(0.1, 1)
)
ggsave(
    filename = file.path(Dir.Exports, "Figure_BoxDEMO_AssocPanels_gg.jpg"),
    plot = AssocPanels_gg2,
    width = 14 * 1.7, height = 16 * 1.5, units = "cm"
)

# ALL-IN-ONE ASSOCIATIONS WITH ENVIRONMENT =================================
FNAME <- file.path(Dir.Concept, "BoxDEMO_AssocsAllinOne_SimResult.RData")

Network_igraph <- Sim.Network(
    n_spec = 6,
    NetworkType = "Association",
    Sparcity = 0,
    MaxStrength = 1, # what absolute value is the maximum link strength
    seed = 42 # seed for randomness
)
E(Network_igraph)$weight <- rep(0, length(E(Network_igraph)))
E(Network_igraph)$weight[1] <- 0.5 # sp1 likes sp2 and vice versa
E(Network_igraph)$weight[length(E(Network_igraph))] <- -0.5 # sp5 hates sp6 and cive versa
## links spanning niches
# E(Network_igraph)$weight[6] <- 1 # sp2 and sp3 like each other
# E(Network_igraph)$weight[11] <- 1 # sp3 and sp5 like each other

CarryingK_vec <- rep(200, length(V(Network_igraph)))
names(CarryingK_vec) <- paste0("Sp_", V(Network_igraph))

Niches_vec <- c(Sp_1 = 5, Sp_2 = 5, Sp_3 = 10, Sp_4 = 10, Sp_5 = 15, Sp_6 = 15)

Effect_Mat <- igraph::as_adjacency_matrix(Network_igraph, attr = "weight")
rownames(Effect_Mat) <- colnames(Effect_Mat) <- names(CarryingK_vec)
Effect_Mat

Initialise_df <- Sim.Initialise(
    n_spec = 6,
    n_individuals = 2e3,
    n_mode = "each",
    Env_range = c(0, 10),
    Trait_means = Niches_vec,
    Trait_sd = 0,
    seed = 42
)
# Plot_IndivsInSpace(Initialise_df)

Env_mat <- Sim.Space(
    x_range = c(0, 20),
    y_range = c(0, 20),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) x,
    y_gradient = function(y) 0
)
# Plot_Environment(Env_mat)

if (file.exists(FNAME)) {
    load(FNAME)
} else {
    SimResult <- Sim.Compute(
        # Demographic parameters
        d0 = d0,
        b0 = b0,
        k_vec = CarryingK_vec,
        ID_df = Initialise_df,

        # Spatial parameters
        env.xy = Env_mat,
        env.sd = 1,
        mig.sd = 0.5,
        mig.top = 0.15,

        # Interaction parameters
        interac.maxdis = 0.5,
        interac.igraph = Network_igraph,
        interac.scale = 1,

        # Simulation parameters
        Sim.t.max = t_max,
        Sim.t.inter = t_inter,
        seed = 42,
        verbose = verbose, # whether to print progress in time as current time
        RunName = "BoxDEMO_Assocs"
    )
    save(SimResult, file = FNAME)
}

last_df <- SimResult[[length(SimResult)]]

p_main <- Plot_Environment(Env_mat) +
    geom_point(
        data = last_df,
        aes(x = X, y = Y, shape = Species),
        fill = "white"
    ) +
    scale_shape_manual(values = 1:length(unique(last_df$Species)))

p_shape <- ggplot(
    last_df,
    aes(x = X, y = Y, shape = Species)
) +
    geom_point() +
    scale_shape_manual(values = 1:length(unique(last_df$Species))) +
    theme_void() +
    theme(legend.position = "right")

shape_legend <- get_legend(p_shape)

p_main_clean <- p_main +
    guides(shape = "none")

SpecinSpace2_gg <- plot_grid(
    p_main_clean,
    shape_legend,
    ncol = 2,
    rel_widths = c(1, 0.1)
)

print(SpecinSpace2_gg)

# INTERACTIONS WITH ENVIRONMENT ============================================
stop("make ring and trophic network")
# maybe make one unaffected primary, one secondary consumer, and one tertiary predator? maybe with weaker links, specify ncihe preference of primary and secondary consumer as increasingly unaligned with environment (do this via negative preferences)

Env_mat <- Sim.Space(
    x_range = c(0, 10),
    y_range = c(0, 10),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) {
        2 * ((x - 5)^2) / 2
    },
    y_gradient = function(y) {
        2 * ((y - 5)^2) / 2
    }
) / 5
Env_mat <- abs(Env_mat - 3)
Plot_Environment(Env_mat)
