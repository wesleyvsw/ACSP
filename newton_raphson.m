%fpt(ind)=5*cos(i)^2-3*sin(3*i)
%derivada = -9*cos(3*i)-10*sen(x)*cos(x)
contador = 0;
x2 = 5;
erro = 1000
while erro > 0.00001 && contador < 100
    x1 =x2
    x2 =  x1 - (5*cos(x1)^2-3*sin(3*x1))/(-9*cos(3*x1)-10*sin(x1)*cos(x1));
    erro = abs(x2-x1)/abs(x2);
    contador = contador+1;
end
x2
contador
