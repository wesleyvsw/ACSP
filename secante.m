%metodo da secante
%fpt(ind)=5*cos(i)^2-3*sin(3*i);
x0 = 0
x1 = 1
contador = 1;
imax = 100;
erro = abs(x1-x0)
erro_max = 0.001
while erro > erro_max
    fx_0 = 5*cos(x0)^2 - 3*sin(3*x0);
    fx_1 = 5*cos(x1)^2 - 3*sin(3*x1);
    x2 = (fx_1*x0 - fx_0*x1)/(fx_1=fx_0);
    x0 = x1;
    x1 = x2;
end
x1
