#!/bin/sh
counter=0
counter2=2
rm -rf ${input}
while [ $counter -lt ${nvars} ]; do
    let counter=counter+1
    let counter2=counter2+1
    head -$counter2 $1 | tail -1 >> ${input}
done
echo ${PID} > solver_script.pid
counter=0
${executable} > executable.log

eval_time=`date +%s`
let exec_time=eval_time-${start_time}
read iter < iteration.number
let iter=iter+1
echo $iter > iteration.number
read -r objective < ${output}
pmod=${printopt}
if [ -f best_objective ]; then
   read bobjective < best_objective
   dif=$(echo "$objective $bobjective" | awk '{if($1 > $2) print 1; else print 2}')
   if [ $dif -eq 2 ]; then
      echo $objective > temp2.txt
      tr -s '\n' ' ' < temp2.txt > temp.txt
      echo >> temp.txt
      cp temp.txt best_objective
      cp ${input} best_solution
      cp iteration.number best_iteration
   fi
else
   echo $objective > temp2.txt
   tr -s '\n' ' ' < temp2.txt > temp.txt   
   echo >> temp.txt
   cp temp.txt best_objective
   cp ${input} best_solution
   cp iteration.number best_iteration
fi
if [ $pmod -gt 0 ]; then
   printmod=$(($iter%$pmod))
   if [ $printmod -eq 0 ]; then
      tr -s '\n' ' ' < iteration.number >> evals.res
      echo $exec_time > temp.txt
      tr -s '\n' ' ' < temp.txt >> evals.res
      tr -s '\n' ' ' < ${input} >> evals.res
      echo $objective > temp.txt
      tr -s '\n' ' ' < temp.txt >> evals.res
      cat best_objective >> evals.res
   fi
fi

echo 1 > $2
echo $objective >> $2

if [ -f mydfopid ]; then
   let endcounter=${limit_evals}+50
   if [ $iter -gt $endcounter ]; then
      kill -9 `cat mydfopid`
   fi
fi

#dif_tol=$(echo "0 ${global_tolerance}" | awk '{if ($1 < $2) print 2; else print 1}')
#if [ $dif_tol -eq 2 ]; then
#   target_solution=`echo "${global_solution}+${global_tolerance}" | bc`
#   dif_obj=$(echo "$objective $target_solution" | awk '{if ($1 > $2) print 1; else print 2}')
#   if [ $dif_obj -eq 2 ]; then
#      kill -9 `cat mydfopid`
#   fi
#fi
