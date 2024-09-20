%Matriz de impedancias
matriz_impedancias = [0   0+0.015i  0   0;0.015i 0 0.0009+0.014i 0; 0 0.0009+0.014i 0 0.015i; 0 0 0.015i 0];
pot_ativa = [4,0,-2,0];
%matriz de gkm
gkm = zeros(4,4);
%parte real da matriz de impedancias
parte_real = real(matriz_impedancias);
%parte imaginaria da matriz de impedâacias
parte_imaginaria = imag(matriz_impedancias);
%Invertendo os elementos da matriz de parte imaginaria
imag_invertida = 1 ./ parte_imaginaria;
% Criando uma máscara lógica para valores infinitos
mascara_infinitos = isinf(imag_invertida);
% Substituindo os valores infinitos por zero
imag_invertida = imag_invertida .* ~mascara_infinitos;
mascara_nan = isnan(imag_invertida);
imag_invertida(mascara_nan) = 0;
%Adicionando a segunda linha em paralelo as contas
imag_invertida(2,3) = 2*imag_invertida(2,3);
imag_invertida(3,2) = 2*imag_invertida(3,2);
%Preenchendo os valores da diagonal principal
[linhas, colunas] = size(imag_invertida);
for i = 1:linhas
    soma_linha = sum(imag_invertida(i, :)) - imag_invertida(i, i);
    imag_invertida(i, i) = soma_linha;
end
%Multiplicando por menos um todos os elementos fora da diagonal principal
for x  = 1:linhas
    for y = 1: colunas
        if x ~= y
            imag_invertida(x,y) = imag_invertida(x,y)*-1;
        end
    end
end
%Metodo linearizado
%Removendo a ultima linha e colunas
imag_invertida = imag_invertida(1:end-1, :);
imag_invertida = imag_invertida(:, 1:end-1);
pot_ativa = pot_ativa(1:end-1); %Removendo o ultimo termo
%angulos theta
angulo = inv(imag_invertida)*pot_ativa';
angulo = [angulo;0];
%Calculando o fluxo de potencia entre as barras
matriz_fluxos = ones(linhas)

for x  = 1:linhas
    for y = 1: colunas
      matriz_fluxos(x,y) = (angulo(x)-angulo(y))/parte_imaginaria(x,y);
    end
end
matriz_fluxos
%Calculando as perdas
for x  = 1:linhas
    for y = 1: colunas
        gkm(x,y) = parte_real(x,y)/(parte_real(x,y)*parte_real(x,y) + parte_imaginaria(x,y)*parte_imaginaria(x,y));
    end
end
gkm(2,3) = gkm(2,3)*2;
gkm(3,2) = gkm(3,2)*2;
perdas23 = gkm(2,3)*(angulo(2)-angulo(3))^2;
%Vetor de perdas nas barras
p_perdas = zeros(4,1);
p_perdas(2) = perdas23/2;
p_perdas(3) = perdas23/2;

p_perdas =p_perdas(1:end-1);
angulo_com_perdas = inv(imag_invertida)*(pot_ativa'-p_perdas);
%Adicionando theta4 = 0 no vetor de angulos
angulo_com_perdas = [angulo_com_perdas;0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%matriz de fluxos com perdas%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matriz_fluxos_com_perdas = ones(linhas)
for x  = 1:linhas
    for y = 1: colunas
      matriz_fluxos_com_perdas(x,y) = (angulo_com_perdas(x)-angulo_com_perdas(y))/parte_imaginaria(x,y);
    end
end
matriz_fluxos_com_perdas
