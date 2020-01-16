function F = hsF37(x)
F(1) = 72 - x(1) - 2*x(2) - 2*x(3);
F(2) = x(1) + 2*x(2) + 2*x(3);
if size(F,1) == 1
  F = F';
end
