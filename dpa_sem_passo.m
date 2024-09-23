%%%%Polos dominantes - resultadso sem controle de passo
delta_s = 1000;
contador = 0;
tau = 1.0;
kg = 10.0;
t1 = 0.005;
s1 = 10i;
s2 = 25i;
s = s1;
gs = kg/(s+t1*s^2+kg*exp(-tau*s))
%derivada = -1*(4000 *exp(s)* exp(-1000 + exp(s)* exp(100 + s))) / (2000 + s*exp(s)*exp(2000 - s))^2
derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2
while abs(delta_s) > 0.01 && contador < 50
  gs = kg/(s+t1*s^2+kg*exp(-tau*s));
  %derivada = -1*(4000 *exp(s)* exp(-1000 + exp(s)* exp(100 + s))) / (2000 + s*exp(s)*exp(2000 - s))^2;
  derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2
  delta_s = gs/derivada
  contador = contador + 1;
  s = s +delta_s;
end
s

