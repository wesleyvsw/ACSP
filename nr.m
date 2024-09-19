matriz_impedancias = [0   0+0.015i  0   0;0.015i 0 0.0009+0.014i 0; 0 0.0009+0.014i 0 0.015i; 0 0 0.015i 0];
matriz_adm_imp = 1./matriz_impedancias;%matriz com o inverso de cada impedancia
b_shunt = 1.2;
% Criando uma máscara lógica para valores infinitos
mascara_infinitos = isinf(matriz_adm_imp);
% Substituindo os valores infinitos por zero
matriz_adm_imp= matriz_adm_imp .* ~mascara_infinitos;
mascara_nan = isnan(matriz_adm_imp);
matriz_adm_imp(mascara_nan) = 0;
%Criando a matriz de admitancia
mat_adm = matriz_adm_imp
mat_adm(2,3) = 2*mat_adm(2,3);
mat_adm(3,2) = 2*mat_adm(3,2)%multiplicando por 2 -->linhas em paralelo
%Adicionando os termos shuntdas barras 2 e 3
mat_adm(2,2) = mat_adm(2,2)+b_shunt;
mat_adm(3,3) = mat_adm(3,3)+b_shunt
%Adicionando os termos da diagonal principal
[linhas, colunas] = size(mat_adm);
for i = 1:linhas
    soma_linha = sum(mat_adm(i, :));
    mat_adm(i, i) = soma_linha;
end
mat_adm
