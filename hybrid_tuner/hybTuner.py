try:
    import ujson as json
except ImportError:
    import json

import numpy as np
import pandas as pd
import argparse
import subprocess
import os
import sys
import time
import random
import shutil
import stat
import math
import matplotlib.pyplot as plt
from mako.template import Template
from hybrid_tuner.dfo_solvers import dfoClass
from hybrid_tuner.pyDirect import pyOpt

class hybClass():
    def __init__(self, args, *pargs, **kwargs):
        np.set_printoptions(precision=3)
        pd.options.display.precision=3
        self.input = 'myin'
        self.output = 'myout'
        self.args = args
        try: 
            os.path.exists(self.args.params)
            self.params = json.load(open(self.args.params))
        except:
            sys.exit("Parameter file does not exist, please provide one.")

        self.nvars = self.params['num_params']
        self.lb = np.asarray(self.params['lower_bounds'])
        self.lb = [float(x) for x in self.lb]
        self.ub = np.asarray(self.params['upper_bounds'])
        self.ub = [float(x) for x in self.ub]
        self.x0 = np.asarray(self.params['starting_point'])
        self.x0 = [float(x) for x in self.x0]

        difL = [x-y for x,y in zip(self.x0, self.lb)]
        for i in range(0, len(difL)):
            if difL[i] < 0:
                self.x0[i] = (self.lb[i] + self.ub[i])/2
        difU = [x-y for x,y in zip(self.x0, self.ub)]
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
        except:
            self.outdir = 'tmp'

        try:
            self.dfo_path = self.params['dfo_path']
        except:
            self.dfo_path = '../..'

        self.printopt = self.params['printopt']
        try:
            self.bandit = self.params['bandit']
        except:
            self.bandit = False

        try:
            self.hybrid = self.params['hybrid']
        except:
            self.hybrid = False

        try:
            self.matlab_path = self.params['matlab_path']
            self.matlab = True
        except:
            self.matlab = False

        self.solver = self.params['solver']

        self.dirpath = os.getcwd()
        self.tdir = self.dirpath + '/' + self.outdir
        self.incumbent = 1000000
        self.s_time = round(time.time())
        self.elapsed = 0
        for i in range(0, len(self.x0)):     
            if self.ints[i] == 1:
                self.x0[i] = round(self.x0[i]) 
        self.time_iter = self.max_cpu
        if(not os.path.exists(self.tdir)):
            os.mkdir(self.tdir)
        if self.bandit:
            err = self.init_bandit()
        if self.hybrid:
            err = self.init_hybrid()
        os.chdir(self.tdir)

    def init_bandit(self):
        try:
            os.path.exists(self.bandit)
            self.bandit_params = json.load(open(self.bandit))
        except:
            sys.exit("Bandit parameter file does not exist, please provide one.")
        self.init = self.bandit_params['init']
        self.frequency = self.bandit_params['frequency']
        self.cParam = self.bandit_params['cParam']
        self.window = self.bandit_params['window']
        self.init_limit = self.bandit_params['init_limit']

        limit = self.frequency
        self.time_iter = self.max_cpu*self.frequency/self.limit_evals
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

        return 0

    def init_hybrid(self):
        try:
            os.path.exists(self.hybrid)
            self.hybrid_params = json.load(open(self.hybrid))
        except:
            sys.exit("Hybrid parameter file does not exist, please provide one.")
        self.hybrid_solver = self.hybrid_params['solver']
        self.init_limit = self.hybrid_params['init_limit']
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

        return 0
        
    def script_setup(self, solver):
        processID = os.getpid()
        template = Template(filename=self.dfo_path + '/solver.script.mako', strict_undefined=True)
        with open(self.tdir + '/' +  str(solver) + '.script', "w") as f:
            f.write(template.render(nvars=self.nvars, PID=processID, executable=self.dirpath + '/' + self.executable, 
                                    start_time=self.s_time, input=self.input, output=self.output, 
                                    printopt=self.printopt, limit_evals=self.limit_evals, 
                                    global_tolerance=self.global_tolerance, global_solution=self.global_solution))
            f.close()
        st = os.stat(self.tdir + '/' + str(solver) + '.script')
        os.chmod(self.tdir + '/' +  str(solver) + '.script', st.st_mode | 0o111)
            
    def matlab_setup(self, solver):
        solver_path = self.dfo_path
        f3 = open(self.tdir + '/matlab_main3.m', 'w')
        f3.write('clear all\n')
        f3.write('quit;\n')
        f3.close()

        var_params = ''
        for i in range(self.nvars):
            var_params += 'bl(' + str(i+1) + ') = ' + str(self.lb[i]) + ';\n' + 'bu(' + str(i+1) + ') = ' + str(self.ub[i]) + ';\n'
            var_params += 'x_0(' + str(i+1) + ') = ' + str(self.x0[i]) + ';\n' + 'ints(' + str(i+1) + ') = ' + str(self.ints[i]) + ';\n'

        template = Template(filename=self.dfo_path + '/matlab_main1.m.mako', strict_undefined=True)
        with open(self.tdir + '/matlab_main1.m', "w") as f:
            f.write(template.render(var_params=var_params, global_solution=self.global_solution, global_tolerance=self.global_tolerance))

        if solver == 'mcs':
            inPair = 'data, x'
        else:
            inPair = 'x, Prob'
        
        fin_loop=''
        for i in range(self.nvars):
            fin_loop += "fprintf(fin, '%30.15f\\n', x(" + str(i+1) + "));\n"

        template = Template(filename=self.dfo_path + '/func_f.m.mako', strict_undefined=True)
        with open(self.tdir + '/func_f.m', "w") as f:
            f.write(template.render(fin_loop=fin_loop, nvars=self.nvars, limit_evals=self.limit_evals, 
                                    executable=self.dirpath+'/'+self.executable, inPair=inPair))
        
        filenames = [self.tdir + '/matlab_main1.m', self.tdir + '/matlab_main2.m', self.tdir + '/matlab_main3.m']
        with open(self.tdir + '/matlab_main.m', 'wb') as wfd:
            for f in filenames:
                with open(f, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
    
    def dfo_iteration(self, solver):
        if solver == 8 or solver == 15:
            self.matlab_setup(solver)
            fun = '-r "addpath(\'' + self.tdir + '\'); matlab_main; exit"'
            cmd = ['matlab', '-nodisplay', '-nosplash', fun]
            with open(self.tdir + '/' + str(solver) +'.results', "w") as outfile:
                code = subprocess.run(cmd, timeout=self.time_iter, stdout=outfile)
        elif solver == 10:
            self.script_setup(solver)
            cmd = [self.dfo_path + '/HOPSPACK_main_serial', 'hopspack_input.in', '&']
            with open(self.tdir + '/' + str(solver) + '.results', "w") as outfile:
                code = subprocess.run(cmd, timeout=self.time_iter, stdout=outfile)
                
    def update_x0(self):
        fn = open('best_objective', 'r')
        curr = float(fn.read())
        fn.close()

        if (curr < self.incumbent):
            self.incumbent = curr
            raw_x0 = np.loadtxt(fname='best_solution', ndmin=1)
            for i in range(0, len(self.x0)):
                if self.ints[i] == 1:
                    self.x0[i] = round(raw_x0[i])

    # The goal of this function is to do a certain number of iterations with a specified solver
    def run_tuning_iteration(self, solver, dfo):
        try:
            dfo.call_dfo(solver)         
        except:
            print('Error producing solver files!')

        self.dfo_iteration(solver)
        fn = open('evals.res', 'r')
        evals = np.genfromtxt(fn)
        try:
            chk = np.size(evals,1) > 1
            df = pd.DataFrame(evals)
        except:
            df = pd.DataFrame(evals).T

        dfI = df.astype({0: 'int32'})
        elapsed = int(dfI.tail(1)[0])            
        dfI[0] = dfI[0] + self.elapsed
        try:
            dfI = dfI.drop([1,5], axis=1) # Drop column related to time and best value for consistency.        
        except:
            dfI = dfI.drop([1], axis=1) # Drop column related to time for consistency.        
        fn.close()
        fsolve = open(str(solver) + '.res', 'a')
        data = dfI.to_string(index=False, header=False) + '\n'
            
        fsolve.write(data)
        fsolve.close()
        fall = open('allEvals.res', 'a')
        fall.write(data)
        fall.close()
        
        self.update_x0()
        try:
            os.remove('evals.res')
        except:
            print('evals.res does not exist')
        return elapsed
        
    def directInit(self):
        nlo = pyOpt(self)
        self.x0, self.incumbent, self.elapsed = nlo.direct(self.init_limit, self.time_iter, self.nvars)
        for i in range(0, len(self.x0)):
            if self.ints[i] == 1:
                self.x0[i] = round(self.x0[i])                            

    def mcsInit(self, dfo):
        dfo.call_dfo(7)
        self.matlab_setup('mcs')
        fun = '-r "addpath(\'' + self.tdir + '\'); matlab_main; exit"'
        cmd = ['matlab', '-nodisplay', '-nosplash', fun]
        with open(self.tdir + '/' + str(solver) +'.results', "w") as outfile:
            code = subprocess.run(cmd, timeout=self.time_iter, stdout=outfile)
        

    def BanditDFO(self, dfo):
        if self.init == 'direct':
            self.directInit()
        elif self.init == 'mcs':
            self.mcsInit(dfo)
        # The following are the list of solvers that we have that accept starting points
        solvers = [8, 10, 15]
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
        
        numer = self.window+self.frequency-2
        limit = numer/self.frequency
        func_evals = np.zeros((len(solvers), int(limit)))
        period = 0
        # The following is the main loop of the bandit function
        while self.elapsed < self.limit_evals:
            fminus = self.run_tuning_iteration(next_solver, dfo) # define this to be an iteration of a specified solver
            self.elapsed = self.elapsed + fminus
            H = self.elapsed
            buff = 0
            if (self.elapsed-self.window > 0):
                buff = self.elapsed - self.window
                H = self.window
            for k in range(len(solvers)):
                idx = solvers[k]
                if next_solver == idx:
                    func_evals[k, round(period % limit)] = fminus
                # Calculate the Ht term:
                counter = 0
                Ht[idx] = sum(func_evals[k])

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
                except:
                    AUC[idx] = 0
                
                if (Ht[idx] != 0):
                    score[idx]=self.cParam*math.sqrt(2*math.log(self.elapsed-buff)/Ht[idx]+AUC[idx])
                else:
                    score[idx] = 100

            #Select next solver
            maxValue = max(score.values())
            keys = [key for key, value in score.items() if value == maxValue]
            next_solver = random.choice(keys)
            period = period + 1            

        #BanditDFO has finished, report results below
        print('BanditDFO has completed!\n')
        print('Best Solution = ' + str(self.incumbent) + ' found after ' + str(self.elapsed) + ' iterations!')       

    def HybridDFO(self, dfo):
        self.directInit()
        elapsed = self.run_tuning_iteration(self.hybrid_solver, dfo) # define this to be an iteration of a specified solver
        self.elapsed = self.elapsed + elapsed

        #HybridDFO has finished, report results below
        print('HybridDFO has completed!\n')
        print('Best Solution = ' + str(self.incumbent) + ' found after ' + str(self.elapsed) + ' iterations!')       

    def SingleSolver(self, dfo):
        solver = int(self.solver)
        f = open(self.tdir + '/iteration.number', 'w')
        f.write('0')
        f.close()

        elapsed = self.run_tuning_iteration(solver, dfo)
        self.elapsed = self.elapsed + elapsed

        #DFO has finished, report results below
        print(str(self.solver) + ' has completed!\n')
        print('Best Solution = ' + str(self.incumbent) + ' found after ' + str(self.elapsed) + ' iterations!')

    
