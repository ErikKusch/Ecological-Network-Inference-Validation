# NO NETWORK SPECIES IN ENVIRONMENT =========================================
Network_igraph <- Sim.Network(
    n_spec = 4,
    NetworkType = "Association",
    Sparcity = 1,
    MaxStrength = 0, # what absolute value is the maximum link strength
    seed = 42 # seed for randomness
)

CarryingK_vec <- c(Sp_1 = 100, Sp_2 = 100, Sp_3 = 100, Sp_4 = 100)

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
    mig.top = 0.25,

    # Interaction parameters
    interac.maxdis = 1,
    interac.igraph = Network_igraph,
    interac.scale = 0,

    # Simulation parameters
    Sim.t.max = t_max,
    Sim.t.inter = t_inter,
    seed = 0,
    verbose = verbose, # whether to print progress in time as current time
    RunName = "BoxDEMO_SpaceOnly"
)

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

final_plot <- plot_grid(
    p_main_clean,
    shape_legend,
    ncol = 2,
    rel_widths = c(1, 0.1)
)

final_plot




# ASSOCIATIONS WITHOUT ENVIRONMENT ========================================



# INTERACTIONS WITHOUT ENVIRONMENT ========================================
