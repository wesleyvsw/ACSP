%valores iniciais
x0 = 0;
x1 = 5;
%iniciando o algoritimo
contador = 1;
erro = 1000;
x = x0;
while erro>0.000001 && contador < 100
    x_anterior = x;
    x = (x0 +x1)/2;
    erro = abs((x-x_anterior)/x);
    fx_0 = 5*cos(x0)^2 - 3*sin(3*x0);
    fx_x = 5*cos(x)^2 - 3*sin(3*x);
    if fx_0*fx_x <0
        x1 = x;
    else
        x0 = x;
    contador = contador +1;
    end
end
x

