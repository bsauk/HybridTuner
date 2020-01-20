
                            SID-PSM Version 1.3
                           =====================

                                 Read me
                               ----------


Welcome to SID-PSM Version 1.3 November, 2014!

This file describes the following topics:


1. System requirements
2. Third party codes
3. Installation
4. Contents of the sid_psm_1.3 directory 
5. The command line
6. The SID-PSM team and web site 
7. Acknowledgments


1. System requirements
----------------------------------------

Your computer should have a Matlab numerical compiler installed. If you plan
to use the current search option implemented in SID-PSM, based on minimum Frobenius
norm  models, then you need a solver for the solution of the trust-region subproblems.
If you chose to use trs function (which is based on dgqt routine from MINPACK2 [1]) 
then you should also have installed a Fortran compiler. Otherwise you need to have 
the Matlab function trust.m, which is part of the Optimization toolbox.

Reference:

[1] J. J. Moré, D. C. Sorensen, K. E. Hillstrom and B. S. Garbow, The MINPACK
    Project, in Sources and Development of Mathematical Software, W. J. Cowell, ed., 
    Prentice-Hall, NJ, 1984, pp. 88--111.


2. Third party codes
----------------------------------------

If you plan to use the search option of SID-PSM then, depending of the choice 
of solver for the trust-region subproblems, one of the following third party 
components is required:

trs_solver = 1: dgqt subroutine from MINPACK2 [1], jointly distributed with
                SID-PSM files.

trs_solver = 0: trust.m Matlab file, from the Optimization toolbox.

Reference:

[1] J. J. Moré, D. C. Sorensen, K. E. Hillstrom and B. S. Garbow, The MINPACK
    Project, in Sources and Development of Mathematical Software, W. J. Cowell, ed., 
    Prentice-Hall, NJ, 1984, pp. 88--111.


3. Installation
----------------------------------------

The code requires approximately 650Kb of hard disk space. If you are running
the code on an Unix/Linux platform, then the following commands should be 
executed to unpack the compressed files:

        gunzip -c sid_psm_1.3.tar.gz | tar xvf -

In a Windows platform use Unzip for the same purpose.

A directory named sid_psm_1.3 will be created, either in Windows or
Unix/Linux operating systems. This directory must be moved to a working
directory in the Matlab tree. Alternatively, the Matlab path could be
updated accordingly.

The SID-PSM parameters can be modified in the file parameters.m
under the sid_psm_1.3 directory. 


4. Contents of the sid_psm_1.3 directory
----------------------------------------

In the directory sid_psm_1.3 there are the following files:

  driver_const.m       A driver for constrained problems.

  driver_unconst.m     A driver for unconstrained problems.

  README.txt           The current file.

  sid_psm_manual.pdf   An user's guide for SID-PSM.

  ---------------------------
  User provided Matlab files:
  ---------------------------

  func_f.m             To compute the value of the objective function.

  func_const.m         To compute the values of the constraint functions.

  grad_const.m         To compute the gradients of the constraint functions.

  -----------------------
  Optimizer Matlab files:
  -----------------------

  domain.m             To check the feasibility of a given point.

  gen.m                To compute the positive spanning set to be used by
                       the pattern search method.

  grad_act.m           To determine the approximated active constraints.

  lambda_poised.m      To compute a \Lambda-poised set from a given
                       set of points.

  match_point.m        To check if a point has been previously evaluated.

  mesh_proc.m          To update the mesh size parameter.

  order_proc.m         To reorder the columns of a given matrix according
                       to the angles between them and a given vector.

  parameters.m         A file with default values for the parameters to be
                       used by SID-PSM.

  proj_ort.m           To compute the orthogonal projection of a given point
                       in a feasible region defined through bounds on the
                       problem variables.

  prune_dir.m          To prune the columns of a given matrix according
                       to the angles between them and a given vector.

  quad_Frob.m          To compute a quadratic interpolation model for the
                       objective function, in the underdetermined case by
                       minimizing the Frobenius norm of the Hessian matrix.

  sid_psm.m            The main program. Applies a pattern search method 
                       to the user's problem.

  simplex_deriv.m      To compute simplex derivatives from a given set of
                       points.

The remaining files in the directory concern the use of the subroutine dgqt 
from MINPACK2 [1].

Reference:

[1] J. J. Moré, D. C. Sorensen, K. E. Hillstrom and B. S. Garbow, The MINPACK
    Project, in Sources and Development of Mathematical Software, W. J. Cowell, ed., 
    Prentice-Hall, NJ, 1984, pp. 88--111.


5. The command line 
----------------------------------------
 
Before running SID-PSM, the user should provide the necessary information
about the objective function in the file func_f.m. Additionally, if 
there are any constraints, the user must supply the corresponding
information in the files func_const.m and grad_const.m. Next, at the
directory sid_psm_1.3, the user should type:

              sid_psm(x_initial,const,output)

where the following data has to be provided:

  - x_initial (the initial point to start the optimizer);

  - const (0 if the problem is unconstrained; 1 if the problem is constrained;
          2 if the constraints are only bounds);

  - output (0-2 variable: 0 - no output; 1 - verbose; 
           2 - very verbose).
 
A description of the algorithmic framework implemented as well as
complementary information about functions, variables and parameters
can be found in the file sid_psm_manual.pdf.

For a quick start, or to confirm that the SID-PSM installation was
completed successfully, the user can run the drivers driver_unconst.m or
driver_const.m.


6. The SID-PSM team and web site 
----------------------------------------

Current SID-PSM team:

   Ana Luisa Custodio (New University of Lisbon).
   Luis Nunes Vicente (University of Coimbra). 

The SID-PSM web site is located at: http://www.mat.uc.pt/sid-psm


7. Acknowledgments
----------------------------------------

SID-PSM team thanks BLAS, LAPACK and MINPACK2 teams for authorizing 
the distribution and use of subroutines jointly with SID-PSM.

