%Variaveis
dt = 0.001;
k = 1.2;
tau = 0.1;
tmax = 100.0;
%Definindo k1 e k2
k1 = (k*dt/2*tau)/(1+(dt/2*tau));
k2 = (1-dt/2*tau)/(1+(dt/2*tau));
% inicializa as variaveis
u(1) = 0;%Vetor de dados de entrada
y(1) = 0;%Vetor de dados de saída
tempo(1) = 0;
contador = 1;
%loop
for i = 0.001:dt:tmax
  if i < 1.0
    u(contador) = 0;
    y(contador) = 0;
  endif
  if i >= 1.0;
    u(contador) = 1;
    y(contador) = k1*u(contador) +k1*u(contador-1) + k2*y(contador-1);
  endif
  tempo(contador) = dt*contador;
  contador = contador + 1;
end
