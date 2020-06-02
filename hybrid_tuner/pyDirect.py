import nlopt
import subprocess


class pyOpt():
    def __init__(self, hybClass):
        self.hybClass = hybClass
        self.iter = 0

    def f(self, x, grad):
        fn = open(self.hybClass.tdir + '/' + self.hybClass.input, "w")
        for i in x:
            fn.write(str(i) + '\n')
        fn.close()
        subprocess.run([self.hybClass.dirpath + '/' +
                        self.hybClass.executable], shell=True)
        fn = open(self.hybClass.tdir + '/' + self.hybClass.output, "r")
        f_out = fn.readline()
        fn.close()
        self.iter = self.iter + 1
        fn = open(self.hybClass.tdir + '/' + 'allEvals.res', "a")
        fn.write(str(self.iter) + '   ' + '   '.join('%.5f' % i for i in x)  +
                 '   ' + '%.5f\n' % float(f_out))
        fn.close()
        return(float(f_out))

    def direct(self, iter_limit, time_limit, nvars):
        opt = nlopt.opt(nlopt.GN_DIRECT_L, nvars)
        x_l = self.hybClass.lb
        x_u = self.hybClass.ub
        x_0 = self.hybClass.x0
        opt.set_lower_bounds(x_l)
        opt.set_upper_bounds(x_u)
        opt.set_maxeval(iter_limit)
        opt.set_maxtime(time_limit)
        opt.set_min_objective(self.f)
        x_opt = opt.optimize(x_0)
        y_opt = opt.last_optimum_value()

        evals = opt.get_numevals()

        return x_opt, y_opt, evals
