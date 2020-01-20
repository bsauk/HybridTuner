global bl
global bu
global global_counter
global global_objective
global global_evaluation
global global_solution
global global_tolerance
global xbest
global start_time

lastn = maxNumCompThreads(1);
start_time = clock;

${var_params}

global_counter = 0;
global_objective = inf;
global_solution = ${global_solution};
global_tolerance = ${global_tolerance};

