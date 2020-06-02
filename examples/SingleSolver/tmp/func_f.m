function f = func_f(x, Prob);

global bl
global bu
global global_counter
global global_objective
global global_evaluation
global global_solution
global global_tolerance
global xbest
global start_time

global_counter = global_counter + 1;
exec_time = etime(clock,start_time);
iten = fopen('iteration.number', 'wt');
fprintf(iten, '%10d\n', global_counter);
fclose(iten);
fin = fopen('myin', 'w');

fprintf(fin, '%30.15f\n', x(1));
fprintf(fin, '%30.15f\n', x(2));
fprintf(fin, '%30.15f\n', x(3));
fprintf(fin, '%30.15f\n', x(4));


fclose(fin);
if ( (x <= bu) & (x >= bl) )
   if (exist('ascale'))
      for i = 1:4
      	  x(i) = x(i)*bscale(i) + ascale(i);
      end;
   end;

   setenv('GFORTRAN_STDIN_UNIT', '5');
   setenv('GFORTRAN_STDOUT_UNIT', '6');
   setenv('GFORTRAN_STDERR_UNIT', '0');
   [status, cmdout] = system('/home/bsauk/Documents/research/5thyear/HybridTuner/examples/SingleSolver/myexec2');
   setenv('GFORTRAN_STDIN_UNIT', '-1') ;
   setenv('GFORTRAN_STDOUT_UNIT', '-1') ;
   setenv('GFORTRAN_STDERR_UNIT', '-1');
   fout = fopen('myout', 'r');
   f = fgetl(fout);
   f = str2num(f);
   fclose(fout);
else
   f = 1e90;
end;

pmod = 1;
if (f < global_objective)
  global_objective = f;
  global_evaluation = global_counter;
  xbest = x;
  gobj = fopen('best_objective', 'wt');
  fprintf(gobj, '%55.8f\n', global_objective);
  fclose(gobj);
  gcou = fopen('best_iteration', 'wt');
  fprintf(gcou, '%10d\n', global_evaluation);
  fclose(gcou);
  gxbest = fopen('best_solution', 'wt');
  fprintf(gxbest, '%10d\n', xbest);
  fclose(gxbest);
end;
if (pmod > 0)
  printmod = mod(global_counter,pmod);
  if (printmod == 0)
    fres = fopen('evals.res', 'a');
    fprintf(fres, '%5d ', global_counter);
    fprintf(fres, '%5d ', exec_time);
    fprintf(fres, '%15.6f ', x);
    fprintf(fres, '%20.8f \n', global_objective);
    fclose(fres);
  end;
end;

if (global_counter > (100))
   exit
end
if (global_tolerance > 0)
   if (global_objective <= (global_solution + global_tolerance) )
      exit
   end
end