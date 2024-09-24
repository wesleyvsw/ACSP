%%%%Polos dominantes - resultadso com controle de passo
delta_s = 1000;
contador_s1 = 0;
contador_s2 = 0;
tau = 1.0;
kg = 10.0;
t1 = 0.005;
s1 = 10i;
s2 = 25i;
s = s1;
derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2;
while abs(delta_s) > 0.001 && contador_s1 < 50
  gs = kg/(s+t1*s^2+kg*exp(-tau*s));
  derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2;
  delta_s = gs/derivada;
  if abs(delta_s) > 0.1;
    delta_s = (delta_s/abs(delta_s)) * 0.1;
end

  s = s + delta_s;
  contador_s1 = contador_s1 + 1;
end
s
contador_s1

s = s2;
delta_s=1000;
derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2;
while abs(delta_s) > 0.001 && contador_s2 < 50
  gs = kg/(s+t1*s^2+kg*exp(-tau*s));
  derivada = -(4000 *exp(s) *(-1000 + exp(s) *(100 + s)))/(2000 + exp(s) *s* (200 + s))^2;
  delta_s = gs/derivada;
  if abs(delta_s) > 0.1;
    delta_s = (delta_s/abs(delta_s)) * 0.1;
end

  s = s + delta_s;
  contador_s2 = contador_s2 + 1;
end
s
contador_s2
