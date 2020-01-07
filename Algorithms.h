class ExpManager;

/**
 * Data transfer functions
 */
void transfer_in(ExpManager* exp_m,bool first_gen);
void transfer_out(ExpManager* exp_m);

void clean(ExpManager* exp_m);
void allocate_next_gen(int nb_indiv);
/**
 * Run all kernel for one generation -- all individuals
 */
void run_a_step_on_GPU(int nb_indiv, double w_max, double selection_pressure, int grid_width, int grid_height, double mutation_rate);
