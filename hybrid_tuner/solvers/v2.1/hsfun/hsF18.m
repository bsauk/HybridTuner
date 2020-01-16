function F = hsF18(x)
F(1) = x(1)*x(2) - 25;
F(2) = x(1)^2 + x(2)^2 - 25;
if size(F,1) == 1
  F = F';
end
