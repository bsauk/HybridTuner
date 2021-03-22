import os


class dfoClass():
    def __init__(self, hybClass):
        self.tune = hybClass 

    def call_dfo(self, idx):
        '''
        Function responsible for calling function to set up DFO scripts

        @type self: Stores hybClass as self.tune
        @param idx: Number of solver scripts to setup
        idx accepts 7, 8, 10, or 15 currently
        @return: Returns the scripts required to run a particular solver
        '''
        if idx == 7:
            chk = self.call_dfo_7()
        elif idx == 8:
            chk = self.call_dfo_8()
        elif idx == 10:
            chk = self.call_dfo_10()
        elif idx == 15:
            chk = self.call_dfo_15()

        return chk

    def call_dfo_7(self):
        '''
        Setup the scripts for solver MCS.
        This requires that self.tune.matlab = True.

        @type self: Stores hybClass as self.tune
        @return: Returns the scripts for MCS
        '''
        f = open(self.tune.tdir + '/matlab_main2.m', 'w')
        f.write("path('" + self.tune.dfo_path + "/mcs', path);\n")
        f.write('n = ' + str(self.tune.nvars) + ';\n')
        f.write("bl = bl';\n")
        f.write("bu = bu';\n")
        f.write('budget = ' + str(self.tune.time_iter) + ';\n')
        f.write("fcn = 'func_f';\n")
        f.write("data = 'func_f';\n")
        f.write('prt = 1;\n')
        f.write('smax = 5*n+10;\n')
        f.write('[x_best, f_best, xmin, fmi, ncall, ncloc] = mcs('
                'fcn, data, bl, bu, prt, smax, budget);\n')
        f.close()

        try:
            os.remove(self.tune.tdir + '/func_f.m')
        except Exception:
            pass

    def call_dfo_8(self):
        '''
        Setup the scripts for solver SID-PSM.
        This requires that self.tune.matlab = True.

        @type self: Stores hybClass as self.tune
        @return: Returns the scripts for SID-PSM
        '''
        f = open(self.tune.tdir + '/matlab_main2.m', 'w')
        f.write("path(path, '" + self.tune.dfo_path + "/sid_psm_1.3');\n")
        f.write('n = ' + str(self.tune.nvars) + ';\n')
        f.write("bl = bl';\n")
        f.write("bu = bu';\n")
        f.write("x_0 = x_0';\n")
        f.write('budget = ' + str(self.tune.time_iter) + ';\n')
        f.write('sid_psm(x_0, 2, 0);\n')
        f.close()
        try:
            os.remove(self.tune.tdir + '/func_f.m')
        except Exception:
            pass

        ffunc = open(self.tune.tdir + '/func_f.m', 'w')
        ffunc.write('function f = func_f(x,Prob);\n')
        ffunc.close()

        fconst = open(self.tune.tdir + '/func_const.m', 'w')
        fconst.write('function [c_const] = func_const(x);\n')
        fconst.write('global bl\n')
        fconst.write('global bu\n')
        fconst.write('c_const = [ ];\n')

        constraint_id = 0
        for i, j in zip(range(0, self.tune.nvars), range(0, 2 *
                                                         self.tune.nvars, 2)):
            fconst.write('c_const(' + str(j+1) + ') = bl(' + str(i+1) + ') - '
                         'x(' + str(i+1) + ');\n')
            fconst.write('c_const(' + str(j+2) + ') = x(' + str(i+1) + ') - '
                         'bu(' + str(i+1) + ');\n')

        fconst.write("c_const = c_const';\n")
        fconst.close()

        fgrad = open(self.tune.tdir + '/grad_const.m', 'w')
        fgrad.write('function [grad_c] = grad_const(x);\n')
        fgrad.write('grad_c = [ ];\n')

        constraint_id = 0
        counter = 0
        while counter < self.tune.nvars:
            counter = counter + 1
            constraint_id = constraint_id + 1
            scounter = 0
            while scounter < self.tune.nvars:
                scounter = scounter + 1
                fgrad.write("grad_c(" + str(scounter) + "," +
                            str(constraint_id) + ") = 0;\n")
            fgrad.write("grad_c(" + str(counter) + "," +
                        str(constraint_id) + ") = -1;\n")
            constraint_id = constraint_id+1
            scounter = 0
            while scounter < self.tune.nvars:
                scounter = scounter + 1
                fgrad.write("grad_c(" + str(scounter) + "," +
                            str(constraint_id) + ") = 0;\n")
            fgrad.write("grad_c(" + str(counter) + "," +
                        str(constraint_id) + ") =  1;\n")
        fgrad.close()

        fp = open(self.tune.tdir + '/parameters.m', 'w')
        fp.write("always = 1;\n")
        fp.write("cache = 0;\n")
        fp.write("economic = 0;\n")
        fp.write("mesh_option = 0;\n")
        fp.write("pss = 2;\n")
        fp.write("regopt = 1;\n")
        fp.write("min_basis = 0;\n")
        fp.write("order_option = 5;\n")
        fp.write("pruning = 0;\n")
        fp.write("search_option = 1;\n")
        fp.write("min_norm = 1;\n")
        fp.write("shessian = 0;\n")
        fp.write("stop_fevals = 1;\n")
        fp.write("stop_alfa = 1;\n")
        fp.write("stop_grad = 0;\n")
        fp.write("stop_iter = 0;\n")
        fp.write("stop_mesh = 1;\n")
        fp.write("store_all = 1;\n")
        fp.write("trs_solver = 0;\n")
        fp.write("alfa = max(1,norm(x_initial,inf));\n")
        fp.write("phi = 1;\n")
        fp.write("theta = 0.5;\n")
        fp.write("fevals_max = " + str(self.tune.frequency) + ";\n")
        fp.write("iter_max = " + str(self.tune.frequency) + ";\n")
        fp.write("tol_alfa = 10^-5;\n")
        fp.write("tol_grad = 10^-5;\n")
        fp.write("epsilon_ini = 10^-4;\n")
        fp.close()

        return True

    def call_dfo_10(self):
        '''
        Setup the scripts for solver HOPSPACK.

        @type self: Stores hybClass as self.tune
        @return: Returns the scripts for HOPSPACK
        '''
        f = open(self.tune.tdir + '/hopspack_input.in', 'w')
        f.write('@ "Problem Definition"\n')
        f.write('"Number Unknowns"'" int " + str(self.tune.nvars) + '\n')
        f.write('"Upper Bounds" vector ' + str(self.tune.nvars) + ' ')
        counter = 0
        while counter < self.tune.nvars:
            f.write(str(self.tune.ub[counter]) + ' ')
            counter = counter + 1
        f.write(' \n')

        f.write('"Lower Bounds" vector ' + str(self.tune.nvars) + ' ')
        counter = 0
        while counter < self.tune.nvars:
            f.write(str(self.tune.lb[counter]) + ' ')
            counter = counter + 1
        f.write(' \n')

        f.write('"Initial X" ')
        f.write("vector " + str(self.tune.nvars) + ' ')
        counter = 0
        for i in range(self.tune.nvars):
            f.write(str(self.tune.x0[i]) + ' ')
        f.write(' \n')
        f.write('"Display"    int 0         \n')
        f.write('@@\n')
        f.write('@ "Evaluator"\n')
        f.write('"Executable Name" string "./10.script"\n')
        f.write('"Input Prefix"   string "hopspackinput"   \n')
        f.write('"Output Prefix"  string "hopspackoutput"  \n')
        f.write('@@\n')
        f.write('@ "Mediator"\n')
        f.write('"Citizen Count" int 1\n')
        f.write('"Maximum Evaluations" int '" " + str(self.tune.frequency) +
                "      \n")
        f.write('"Display"   int 0      \n')
        f.write('@@\n')
        f.write('@ "Citizen 1"\n')
        f.write('"Type" string "GSS"\n')
        f.write('"Step Tolerance" double 0.001\n')
        f.write('@@\n')

        fcall = open(self.tune.tdir + '/10.call', 'w')
        fcall.write(self.tune.dfo_path +
                    '/HOPSPACK_main_serial hopspack_input.in > 10.results &')
        fcall.close()
        f.close()

        return True

    def call_dfo_15(self):
        '''
        Setup the scripts for solver SNOBFIT.
        Requires that self.tune.matlab = True.

        @type self: Stores hybClass as self.tune
        @return: Returns the scripts for SNOBFIT
        '''
        f = open(self.tune.tdir + '/matlab_main2.m', 'w')
        f.write("path(path, '" + self.tune.dfo_path + "/v2.1');\n")
        f.write("file = 'test';\n")
        f.write("fcn = 'func_f';\n")
        f.write('fac = 0;\n')
        f.write('ncall = ' + str(self.tune.frequency) + ';\n')
        f.write('u = bl;\n')
        f.write('v = bu;\n')
        f.write('n = ' + str(self.tune.nvars) + ';\n')
        f.write('npoint = 1;\n')
        f.write('nreq = n + 6;\n')
        f.write('x = x_0;\n')
        f.write("dx = (v-u)'*1.e-5;\n")
        f.write("for i=1:n\n")
        f.write("   if(dx(i) > 1.e-5)\n")
        f.write("      dx(i) = 1.e-5\n")
        f.write("   end\n")
        f.write("end\n")
        f.write("p = 0.5;\n")
        f.write("prt = 0;\n")

        f.write('for j=1:npoint\n')
        f.write('   f(j,:) = [feval(fcn, x(j,:))+fac*randn max(sqrt(eps), '
                '3*fac)];\n')
        f.write('end\n')
        f.write('ncall0 = npoint;\n')
        f.write("params = struct('bounds', {u,v}, 'nreq', nreq, 'p', p);\n")
        f.write('% repeated calls to Snobfit\n')
        f.write('while ncall0 < ncall %repeat till ncall function values '
                'are reached\n')
        f.write('   if ncall0 == npoint % intial call\n')
        f.write('      [request, xbest, fbest] = '
                'snobfit(file, x, f, params, dx);\n')
        f.write('      ncall0, xbest, fbest\n')
        f.write('   else\n')
        f.write('      [request, xbest, fbest] = '
                'snobfit(file, x, f, params);\n')
        f.write('   end\n')
        f.write('   if prt>0, request, end\n')
        f.write('   clear x\n')
        f.write('   clear f\n')
        f.write('   for j=1:size(request,1)\n')
        f.write('      x(j,:) = request(j,1:n);\n')
        f.write('      f(j,:) = [feval(fcn, x(j,:))+fac*randn '
                'max(sqrt(eps), 3*fac)];\n')
        f.write('   end\n')
        f.write('   ncall0 = ncall0 + size(f,1);\n')
        f.write('   [fbestn jbest] = min(f(:,1));\n')
        f.write('   if fbestn < fbest\n')
        f.write('      fbest = fbestn\n;')
        f.write('      xbest = x(jbest, :);\n')
        f.write('      ncall0, xbest, fbest % Display Current number of '
                'function values, best point and function values if '
                'fbest has changed\n')
        f.write('   end\n')
        f.write('end\n')
        f.write('ncall0, xbest, fbest %show number of function values, '
                'best point and function value\n')

        try:
            os.remove(self.tune.tdir + '/func_f.m')
        except Exception:
            pass

        ffunc = open(self.tune.tdir + '/func_f.m', 'w')
        ffunc.write('function f = func_f(x,Prob);\n')
        ffunc.close()
        f.close()

        return True
