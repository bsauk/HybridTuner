function [c_const] = func_const(x);
global bl
global bu
c_const = [ ];
c_const(1) = bl(1) - x(1);
c_const(2) = x(1) - bu(1);
c_const(3) = bl(2) - x(2);
c_const(4) = x(2) - bu(2);
c_const(5) = bl(3) - x(3);
c_const(6) = x(3) - bu(3);
c_const(7) = bl(4) - x(4);
c_const(8) = x(4) - bu(4);
c_const = c_const';
