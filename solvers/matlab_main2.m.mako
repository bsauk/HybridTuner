path(path,'${path_tomlab}');
startup  %include TOMLAB directories to path
n = ${nvars}
Name = 'noname';
x_opt = [];
x_min = bl';
x_max = bu';
bl = transpose(bl);
bu = transpose(bu);
x_0 = transpose(x_0);
Prob = glcAssign('func_f', bl, bu, Name, [], [], [], ...
                 [], [], [], x_0, ...
                 [], [], [], [], ...
                  [], [], [], [], []);
Prob.optParam.MaxIter = ${limit_evals};
Prob.optParam.MaxFunc = ${limit_evals}
Prob.MaxCPU = ${max_cpu};
Prob.LGO.options.G_maxfct = ${limit_evals};
Prob.GO.MaxFunc = ${limit_evals};
Prob.CGO.rbfType = 2;
Prob.x_0 = x_0;
Prob.MIP.IntVars = ints;
Result1 = tomRun('glcDirect', Prob, 1);