function [c_const] = func_const(x);
global bl
global bu
c_const = [ ];
c_const(1) = bl(1) - x(1);
c_const(2) = x(1) - bu(1);
c_const(3) = bl(2) - x(2);
c_const(4) = x(2) - bu(2);
c_const = c_const';
