// General A/B cores
X(  a_core_diameter,                md::scalar,     0.2                     )
X(  b_core_diameter,                md::scalar,     0.2                     )
X(  a_core_repulsion,               md::scalar,     2.0                     )
X(  b_core_repulsion,               md::scalar,     2.0                     )

// Chromatin beads
X(  chromatin_bond_spring,          md::scalar,     0.1                     )
X(  chromatin_bond_length,          md::scalar,     0.2                     )
X(  chromatin_mobility,             md::scalar,     1.0                     )

X(  a_core_bond_spring,             md::scalar,     0.0                     )
X(  a_core_bond_length,             md::scalar,     0.0                     )
X(  b_core_bond_spring,             md::scalar,     0.0                     )
X(  b_core_bond_length,             md::scalar,     0.0                     )

X(  a_core_2nd_bond_spring,         md::scalar,     0.0                     )
X(  b_core_2nd_bond_spring,         md::scalar,     0.0                     )

// Nucleolar particles
X(  nucleolus_sidebeads,            md::index,      2                       )
X(  nucleolus_a_factor,             md::scalar,     5                       )
X(  nucleolus_b_factor,             md::scalar,     5                       )
X(  nucleolus_bond_spring,          md::scalar,     5.0                     )
X(  nucleolus_bond_length,          md::scalar,     0.0                     )
X(  nucleolus_droplet_energy,       md::scalar,     0.0                     )
X(  nucleolus_droplet_decay,        md::scalar,     0.2                     )
X(  nucleolus_droplet_cutoff,       md::scalar,     0.4                     )
X(  nucleolus_mobility,             md::scalar,     1.0                     )

// Nuclear membrane
X(  wall_init_semiaxes,             md::vector,     V(1.0, 1.0, 1.0)        )
X(  wall_semiaxes_spring,           md::vector,     V(1.0e4, 1.0e4, 1.0e4)  )
X(  wall_packing_spring,            md::scalar,     5000                    )
X(  wall_a_factor,                  md::scalar,     5                       )
X(  wall_b_factor,                  md::scalar,     5                       )
X(  wall_mobility,                  md::scalar,     1.0e-4                  )

// Bead scaling parameters
X(  bead_scale_init,                md::scalar,     1.0                     )
X(  bead_scale_tau,                 md::scalar,     1.0                     )
X(  bond_scale_init,                md::scalar,     1.0                     )
X(  bond_scale_tau,                 md::scalar,     1.0                     )

// Initialization stage (spindle and packing)
X(  init_coarse_graining,           md::index,      100                     )
X(  init_bead_diameter,             md::scalar,     0.2                     )
X(  init_bead_repulsion,            md::scalar,     5.0                     )
X(  init_bond_length,               md::scalar,     0.2                     )
X(  init_bond_spring,               md::scalar,     500.0                   )
X(  init_bend_energy,               md::scalar,     0.0                     )
X(  init_spindle_spring,            md::scalar,     1.0                     )
X(  init_spindle_point,             md::point,      V(0, 0, 0)              )
X(  init_packing_radius,            md::scalar,     1.0                     )
X(  init_packing_spring,            md::scalar,     0.0                     )
X(  init_start_point,               md::point,      V(5, 0, 0)              )
X(  init_start_stddev,              md::scalar,     1.0                     )
X(  init_mobility,                  md::scalar,     1.0                     )
X(  init_temperature,               md::scalar,     0.1                     )
X(  init_timestep,                  md::scalar,     1e-4                    )
X(  init_spacestep,                 md::scalar,     0                       )
X(  init_spindle_steps,             md::step,       10000                   )
X(  init_packing_steps,             md::step,       10000                   )
X(  init_sampling_interval,         md::step,       1000                    )
X(  init_logging_interval,          md::step,       1000                    )
X(  init_refinement_method,         std::string,    "spline"                )

// Relaxation stage
X(  relaxation_temperature,         md::scalar,     1.0                     )
X(  relaxation_timestep,            md::scalar,     1.0e-5                  )
X(  relaxation_spacestep,           md::scalar,     0                       )
X(  relaxation_steps,               md::step,       10000                   )
X(  relaxation_sampling_interval,   md::step,       100                     )
X(  relaxation_logging_interval,    md::step,       100                     )

// Interphase stage
X(  interphase_temperature,         md::scalar,     1.0                     )
X(  interphase_timestep,            md::scalar,     1e-5                    )
X(  interphase_spacestep,           md::scalar,     0                       )
X(  interphase_steps,               md::step,       10000                   )
X(  interphase_sampling_interval,   md::step,       1000                    )
X(  interphase_logging_interval,    md::step,       100                     )

X(  contactmap_distance,            md::scalar,     0.4                     )
X(  contactmap_update_interval,     md::step,       100                     )
X(  contactmap_thinning_rate,       md::step,       100                     )

// Random seed values for spindle initialization and relaxation/interphase stages
X(  spindle_seed,                   std::uint64_t,  0                       )
X(  interphase_seed,                std::uint64_t,  0                       )
