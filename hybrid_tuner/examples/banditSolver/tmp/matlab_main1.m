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

bl(1) = 3.0;
bu(1) = 10.0;
x_0(1) = 3.0;
ints(1) = 1;
bl(2) = 3.0;
bu(2) = 10.0;
x_0(2) = 3.0;
ints(2) = 1;


global_counter = 0;
global_objective = inf;
global_solution = 0;
global_tolerance = 0.001;

