A=[-1;1];
b=[2;-4];
cvx_begin
variable x(1)
dual variable y
minimize (x*x+1);
A*x+b <= 0 : y;
cvx_end