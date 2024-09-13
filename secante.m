%metodo da secante
%fpt(ind)=5*cos(i)^2-3*sin(3*i);
%chutes iniciais
x0 = 0;
x1 = 5;
%inicio da funcao
contador = 1;
imax = 100;
erro = 1000;
erro_max = 0.0001;
xl1 = x0;
x_calc = x1;
while erro > erro_max  && contador < 100
    xl2 = xl1;
    xl1 = x_calc;
    fx_l2 = 5*cos(xl1)^2 - 3*sin(3*xl2);
    fx_l1 = 5*cos(xl1)^2 - 3*sin(3*xl1);
    x_calc = (fx_l1*(xl2-xl1))/(fx_l2-fx_l1);
    erro = abs(x_calc-xl1)/x_calc;
    contador = contador +1;
end
x_calc
