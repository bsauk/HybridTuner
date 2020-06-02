try:
    import ujson as json
except ImportError:
    import json

import numpy as np
import pandas as pd
import subprocess
import os
import sys
import time
import random
import shutil
import math
import matplotlib.pyplot as plt
from mako.template import Template
from hybrid_tuner.pyDirect import pyOpt

class hybClass():
    '''
    This is the definition of the hybClass module.
    This module contains Bandit DFO and Hybrid DFO.
    The only requirement for this module is a defined myparams.json file.
    The params file should provide:
    a. num_params: Number of hyper parameters
    b. lower_bounds: Variable lower bounds provided as an array with length num_params
    c. upper_bounds: Variable upper bounds provided as an array with length num_params
    d. starting_point: Variable starting points provided as an array with length num_params
    e. var_type: Variable type, 0 is cont. and 1 is integer provided as array with length num_params
    f. max_iterations: Number of iterations to the black-box function
    g. global_tol: Stops within tolerance of the optimal solution if provided
    h. cpu_limit: Max amount of time to spend on the search provided in seconds
    i. global_target: Optimal solution if one is known
    j. executable: Name of executable file or black-box to call, see myexec
    k. outdir: Location of output directory, defaults to ./tmp
    l. bandit (optional): Location of bparams.json, file of bandit params
    m. hybrid (optional): Location of hparams.json, file of hybrid params
    n. printopt: 0 or 1, if 1 will print intermediate values to evals.res
    o. solver: DFO solver to use initially, accepted options are 8, 10, 15
    '''
    def __init__(self, args, *pargs, **kwargs):
        np.set_printoptions(precision=3)
        pd.options.display.precision = 3
        self.input = 'myin'
        self.output = 'myout'
        self.args = args
        if os.path.exists(self.args.params):
            self.params = json.load(open(self.args.params))
        else:
            sys.exit("Parameter file does not exist, please provide one.")

        self.nvars = self.params['num_params']
        self.lb = np.asarray(self.params['lower_bounds'])
        self.lb = [float(x) for x in self.lb]
        self.ub = np.asarray(self.params['upper_bounds'])
        self.ub = [float(x) for x in self.ub]
        self.x0 = np.asarray(self.params['starting_point'])
        self.x0 = [float(x) for x in self.x0]

        difL = [x-y for x, y in zip(self.x0, self.lb)]
        for i in range(0, len(difL)):
            if difL[i] < 0:
                self.x0[i] = (self.lb[i] + self.ub[i])/2
        difU = [x-y for x, y in zip(self.x0, self.ub)]
        for i in range(0, len(difU)):
            if difU[i] > 0:
                self.x0[i] = (self.ub[i] + self.lb[i])/2

        self.ints = self.params['var_type']
        self.frequency = self.params['max_iterations']
        self.limit_evals = self.params['max_iterations']
        self.window = self.params['max_iterations']
        self.global_tolerance = self.params['global_tol']
        self.cParam = 0.05
        self.init = 0
        self.init_limit = 0

        self.max_cpu = self.params['cpu_limit']
        self.global_solution = self.params['global_target']
        self.executable = self.params['executable']
        try:
            self.outdir = self.params['outdir']
        except KeyError:
            self.outdir = 'tmp'

        self.dfo_path = os.path.abspath(os.path.join(
            os.path.dirname(__file__), "../solvers"))

        self.printopt = self.params['printopt']
        try:
            self.bandit = self.params['bandit']
        except KeyError:
            self.bandit = False

        try:
            self.hybrid = self.params['hybrid']
        except KeyError:
            self.hybrid = False

        self.matlab_path = shutil.which('matlab')
        if self.matlab_path:
            self.matlab = True
        else:
            self.matlab = False

        self.solver = self.params['solver']

        self.dirpath = os.getcwd()
        self.tdir = self.dirpath + '/' + self.outdir
        self.incumbent = 1000000.00
        self.s_time = round(time.time())
        self.elapsed = 0
        for i in range(0, len(self.x0)):
            if self.ints[i] == 1:
                self.x0[i] = round(self.x0[i])
        self.time_iter = self.max_cpu
        if(not os.path.exists(self.tdir)):
            os.mkdir(self.tdir)
        if self.bandit:
            self.init_bandit()
        if self.hybrid:
            self.init_hybrid()
        os.chdir(self.tdir)

    def init_bandit(self):
        '''
        Initializes the bandit parameters for Bandit DFO

        @type self.bandit: Location of bandit parameter file
        @return: Values stored in self for Bandit DFO
        '''
        try:
            os.path.exists(self.bandit)
            self.bandit_params = json.load(open(self.bandit))
        except KeyError:
            sys.exit("Bandit file does not exist, please provide one.")
        self.init = self.bandit_params['init']
        self.frequency = self.bandit_params['frequency']
        self.cParam = self.bandit_params['cParam']
        self.window = self.bandit_params['window']
        self.init_limit = self.bandit_params['init_limit']
        self.time_iter = self.max_cpu*self.frequency/self.limit_evals
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

    def init_hybrid(self):
        '''
        Initializes the bandit parameters for Hybrid DFO

        @type self.hybrid: Location of hybrid parameter file
        @return: Values stored in self for Hybrid DFO
        '''
        try:
            os.path.exists(self.hybrid)
            self.hybrid_params = json.load(open(self.hybrid))
        except KeyError:
            sys.exit("Hybrid file does not exist, please provide one.")
        self.hybrid_solver = self.hybrid_params['solver']
        self.init_limit = self.hybrid_params['init_limit']
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

    def script_setup(self, solver):
        '''
        Initializes the DFO solvers that execute a script

        @type self: Parameters defined in __init__
        @param solver: Currently only accepts solver 10
        @return: Setup the solver file for execution
        '''
        processID = os.getpid()
        template = Template(filename=self.dfo_path + '/solver.script.mako',
                            strict_undefined=True)
        with open(self.tdir + '/' + str(solver) + '.script', "w") as f:
            f.write(template.render(nvars=self.nvars, PID=processID,
                                    executable=self.dirpath + '/' +
                                    self.executable, start_time=self.s_time,
                                    input=self.input, output=self.output,
                                    printopt=self.printopt,
                                    limit_evals=self.frequency,
                                    global_tolerance=self.global_tolerance,
                                    global_solution=self.global_solution))
            f.close()
        st = os.stat(self.tdir + '/' + str(solver) + '.script')
        os.chmod(self.tdir + '/' + str(solver) + '.script', st.st_mode | 0o111)

    def matlab_setup(self, solver):
        '''
        Initialize DFO solvers that use MATLAB,
        only works if self.matlab = True

        @type self: Parameters defined in __init__
        @param solver: Accepts number for MATLAB solvers (8 or 15)
        @return: Matlab scripts to run solver 8 or 15
        '''
        f3 = open(self.tdir + '/matlab_main3.m', 'w')
        f3.write('clear all\n')
        f3.write('quit;\n')
        f3.close()

        var_params = ''
        for i in range(self.nvars):
            lb = str(self.lb[i])
            ub = str(self.ub[i])
            x0 = str(self.x0[i])
            ints = str(self.ints[i])
            idx = str(i+1)
            var_params = (var_params + 'bl(' + idx + ') = ' + lb + ';\n' +
                          'bu(' + idx + ') = ' + ub + ';\n')
            var_params = (var_params + 'x_0(' + idx + ') = ' + x0 + ';\n' +
                          'ints(' + idx + ') = ' + ints + ';\n')

        template = Template(filename=self.dfo_path + '/matlab_main1.m.mako',
                            strict_undefined=True)
        with open(self.tdir + '/matlab_main1.m', "w") as f:
            f.write(template.render(var_params=var_params,
                                    global_solution=self.global_solution,
                                    global_tolerance=self.global_tolerance))

        if solver == 7:
            inPair = 'data, x'
        else:
            inPair = 'x, Prob'
        fin_loop = ''
        for i in range(self.nvars):
            fin_loop += "fprintf(fin, '%30.15f\\n', x(" + str(i+1) + "));\n"

        template = Template(filename=self.dfo_path + '/func_f.m.mako',
                            strict_undefined=True)
        with open(self.tdir + '/func_f.m', "w") as f:
            f.write(template.render(fin_loop=fin_loop,
                                    nvars=self.nvars,
                                    limit_evals=self.frequency,
                                    executable=self.dirpath+'/'
                                    + self.executable, inPair=inPair))
            filenames = [self.tdir + '/matlab_main1.m',
                         self.tdir + '/matlab_main2.m',
                         self.tdir + '/matlab_main3.m']
        with open(self.tdir + '/matlab_main.m', 'wb') as wfd:
            for f in filenames:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)

    def dfo_iteration(self, solver):
        '''
        Executes one call to DFO solver 7, 8, 10, or 15

        @type self: Parameters setup in __init__
        @param solver: Number of DFO solver to call
        @return: Performs one call to the DFO solver writes to
        evals.res if printopt = 1
        '''
        if solver == 7 or solver == 8 or solver == 15:
            self.matlab_setup(solver)
            fun = '-r "addpath(\'' + self.tdir + '\'); matlab_main; exit"'
            cmd = ['matlab', '-nodisplay', '-nosplash', fun]
            with open(self.tdir + '/' + str(solver)
                      + '.results', "w") as outfile:
                subprocess.run(cmd, timeout=self.time_iter, stdout=outfile)
        elif solver == 10:
            self.script_setup(solver)
            cmd = [self.dfo_path + '/HOPSPACK_main_serial',
                   'hopspack_input.in', '&']
            with open(self.tdir + '/' + str(solver)
                      + '.results', "w") as outfile:
                subprocess.run(cmd, timeout=self.time_iter, stdout=outfile)

    def update_x0(self):
        '''
        Updates the best solution and parameters to obtain it
        after running a DFO iteration

        @type self: Parameters setup in __init__
        @return: Update to self.x0 and self.incumbent
        '''
        fn = open('best_objective', 'r')
        curr = float(fn.read())
        fn.close()

        if (curr < self.incumbent):
            self.incumbent = curr
            raw_x0 = np.loadtxt(fname='best_solution', ndmin=1)
            for i in range(0, len(self.x0)):
                if self.ints[i] == 1:
                    self.x0[i] = round(raw_x0[i])

    def run_tuning_iteration(self, solver, dfo):
        '''
        Master function to call one iteration of tuning.

        @type self: Parameters setup in __init__
        @param solver: Number of DFO solver to invoke
        @param dfo: dfoClass see dfo_solvers for more info.
        @return: Updates to self.x0, self.incumbent
        and history of the iterations if printopt = 1
        '''
        try:
            dfo.call_dfo(solver)
        except KeyError:
            print('Error producing solver files!')

        self.dfo_iteration(solver)
        fn = open('evals.res', 'r')
        evals = np.genfromtxt(fn)
        if np.size(evals, 1) > 1:
            df = pd.DataFrame(evals)
        else:
            df = pd.DataFrame(evals).T

        dfI = df.astype({0: 'int32'})
        elapsed = int(dfI.tail(1)[0])
        dfI[0] = dfI[0] + self.elapsed
        if (len(dfI.columns) > self.nvars+3):
            dfI = dfI.drop([self.nvars+3], axis=1) # Drop column with best iteration and only keep current iteration
        dfI = dfI.drop([1], axis=1)  # Drop odd column
        fn.close()
        fsolve = open(str(solver) + '.res', 'a')
        data = dfI.to_string(index=False, header=False) + '\n'
        fsolve.write(data)
        fsolve.close()
        fall = open('allEvals.res', 'a')
        fall.write(data)
        fall.close()
        self.update_x0()
        self.elapsed += elapsed
        try:
            os.remove('evals.res')
        except KeyError:
            print('evals.res does not exist')

        return elapsed

    def hybInit(self, dfo):
        '''
        Execute the global DFO strategy to init Bandit DFO and Hybrid DFO.
        Only works for DIRECT (from NLOPT) or MCS if self.matlab = True.

        @type self: Parameters setup in __init__
        @type dfo: dfoClass that creates DFO scripts
        @return: Executes the initialization of Bandit DFO and Hybrid DFO
        '''
        if self.init == 'direct':
            nlo = pyOpt(self)
            self.x0, self.incumbent, self.elapsed = nlo.direct(self.init_limit,
                                                               self.time_iter,
                                                               self.nvars)
            for i in range(0, len(self.x0)):
                if self.ints[i] == 1:
                    self.x0[i] = round(self.x0[i])
        elif self.init == 'mcs':
            if self.matlab:
                self.run_tuning_iteration(7, dfo)
            else:
                print('MCS cannot be used if MATLAB is not installed!')

    def BanditDFO(self, dfo):
        '''
        Performs Bandit DFO. This will run different solvers until terminating.

        @type self: Parameters setup in __init__
        @type dfo: dfoClass that creates DFO scripts
        @return: Performs Bandit DFO and outputs results interpretable
        by self.visualizeResults
        '''
        self.hybInit(dfo)
        if self.matlab:
            # The following are the list of solvers that we that accept x0
            solvers = [8, 10, 15]
        else:
            solvers = [10]  # Bandit = Hybrid if Matlab is not installed
        if self.solver:
            next_solver = int(self.solver)
        else:
            next_solver = solvers[random.randint(0, len(solvers)-1)]
        AUC = {}
        score = {}
        Ht = {}
        for i in solvers:
            AUC[i] = 2
            score[i] = 2
            Ht[i] = 0

        func_evals = pd.DataFrame([0, 0, 0])
        period = 0
        # The following is the main loop of the bandit function
        while self.elapsed < self.limit_evals:
            fminus = self.run_tuning_iteration(next_solver, dfo)
            if self.incumbent == self.global_solution:
                break
            H = self.elapsed
            buff = 0
            if (self.elapsed-self.window > 0):
                buff = self.elapsed - self.window
                H = self.window
            for k in range(len(solvers)):
                idx = solvers[k]
                if next_solver == idx:
                    func_evals.loc[k, round(period)] = fminus
                # Calculate the Ht term:
                Ht[idx] = sum(func_evals.loc[k])
                # Append a new column to the dataframe if another iteration
                period += 1
                func_evals[period] = 0

                # Now we calculate Vti
                try:
                    fn = open(str(idx) + '.res', 'r')
                    data = np.genfromtxt(fn)
                    small = data[buff][-1]
                    Vti = []
                    for j in range(buff, len(data)):
                        if (data[j][-1] < small):
                            Vti.append(j-buff)
                    lhs = 2/(H*(H+1))
                    rhs = sum(Vti)
                    AUC[idx] = lhs*rhs
                except FileNotFoundError:
                    AUC[idx] = 0

                if (Ht[idx] != 0):
                    score[idx] = self.cParam*math.sqrt(2*math.log(
                        self.elapsed-buff)/Ht[idx]+AUC[idx])
                else:
                    score[idx] = 100

            # Select next solver
            maxValue = max(score.values())
            keys = [key for key, value in score.items() if value == maxValue]
            next_solver = random.choice(keys)

    def HybridDFO(self, dfo):
        '''
        Performs Hybrid DFO. This will run different solvers until terminating.

        @type self: Parameters setup in __init__
        @type dfo: dfoClass that creates DFO scripts
        @return: Performs Hybrid DFO and outputs results interpretable
        by self.visualizeResults
        '''
        self.hybInit(dfo)
        self.run_tuning_iteration(self.hybrid_solver, dfo)

    def SingleSolver(self, dfo):
        '''
        Performs DFO with one solver. This will run until terminating.

        @type self: Parameters setup in __init__
        @type dfo: dfoClass that creates DFO scripts
        @return: Performs DFO and outputs results interpretable
        by self.visualizeResults
        '''
        solver = int(self.solver)
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

        self.run_tuning_iteration(solver, dfo)

    def visualizeResults(self, method):
        '''
        Produces output plots and summarizes results from tuning experiments.

        @type self: Parameters setup in __init__ and generated from DFO
        @type method: Name of method that was invoked
        accepts 'Bandit', 'Hybrid', or 'Single'
        @return: Produces .jpg to display tuning performance over time
        and summarizes results as output
        '''
        if method == 'Bandit':
            print('BanditDFO has completed!')
            print('Best Solution = ' + str(self.incumbent) + ' found after ' +
                  str(self.elapsed) + ' iterations!')
        elif method == 'Hybrid':
            print('HybridDFO has completed!')
            print('Best Solution = ' + '%.5f' % self.incumbent + ' found after ' +
                  str(self.elapsed) + ' iterations!')
#            print('Best Solution = ' + str(self.incumbent) + ' found after ' +
#                  str(self.elapsed) + ' iterations!')
        elif method == 'Single':
            print(str(self.solver) + ' has completed!')
            print('Best Solution = ' + str(self.incumbent) + ' found after ' +
                  str(self.elapsed) + ' iterations!')

        fn = open('allEvals.res', 'r')
        data = np.genfromtxt(fn)
        xtmp = data[:, 0]
        ytmp = data[:, -1]
        y0 = ytmp[0]
        y = []
        x = xtmp.tolist()
        for i in ytmp:
            if i <= y0:
                y0 = i
            y.append(y0)

        plt.scatter(x, y)
        plt.xlabel('Iteration Number')
        plt.ylabel('Objective Value')
        plt.savefig('evalsPlot')
