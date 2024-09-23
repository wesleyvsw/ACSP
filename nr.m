matriz_impedancias = [0   0+0.015i  0   0;0.015i 0 0.0009+0.014i 0; 0 0.0009+0.014i 0 0.015i; 0 0 0.015i 0];
num_barras = 4;
matriz_adm_imp = 1./matriz_impedancias;%matriz com o inverso de cada impedancia
b_shunt = 1.2i;
vetor tensao = [1;1;1;1];
vetor_angulo = [0;0;0;0];
vetor_v_theta = [vetor_angulo(1);vetor_angulo(2);vetor_angulo(3);vetor_tensão(3)];%[theta1,theta2,theta3,v3]
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
mat_adm(3,3) = mat_adm(3,3)+b_shunt;
%Adicionando os termos da diagonal principal
[linhas, colunas] = size(mat_adm);
for i = 1:linhas
    soma_linha = sum(mat_adm(i, :));
    mat_adm(i, i) = soma_linha;
end
%Multiplicando por menos um todos os elementos fora da diagonal principal
for x  = 1:linhas
    for y = 1: colunas
        if x ~= y
            mat_adm(x,y) = mat_adm(x,y)*-1;
        end
    end
end
%parte real da matriz de admitancias
parte_real = real(matriz_adm);
%parte imaginaria da matriz de admitancias
parte_imaginaria = imag(mat_adm);
%%%%%%%%%%%%%%%%%%%%%%%%METODO DE NEWTON%%%%%%%%%%%%%%%%%%%%%%
%POTENCIAS Ativas
p1_calc =  vetor_tensao(1)*vetor_tensao(1)*(parte_real(1,1)*cos(vetor_angulo(1) - vetor_angulo(1))+ parte_imaginaria(1,1)*sin(vetor_angulo(1) - vetor_angulo(1)))+vetor_tensao(1)*vetor_tensao(2)*(parte_real(1,2)*cos(vetor_angulo(1) - vetor_angulo(2))+ parte_imaginaria(1,2)*sin(vetor_angulo(1) - vetor_angulo(2)))+vetor_tensao(1)*vetor_tensao(3)*(parte_real(1,3)*cos(vetor_angulo(1) - vetor_angulo(3))+ parte_imaginaria(1,3)*sin(vetor_angulo(1) - vetor_angulo(3)))+vetor_tensao(1)*vetor_tensao(4)*(parte_real(1,4)*cos(vetor_angulo(1) - vetor_angulo(4))+ parte_imaginaria(1,4)*sin(vetor_angulo(1) - vetor_angulo(4)));
p2_calc =  vetor_tensao(2)*vetor_tensao(1)*(parte_real(2,1)*cos(vetor_angulo(2) - vetor_angulo(1))+ parte_imaginaria(2,1)*sin(vetor_angulo(2) - vetor_angulo(1)))+vetor_tensao(2)*vetor_tensao(2)*(parte_real(2,2)*cos(vetor_angulo(2) - vetor_angulo(2))+ parte_imaginaria(2,2)*sin(vetor_angulo(2) - vetor_angulo(2)))+vetor_tensao(2)*vetor_tensao(3)*(parte_real(2,3)*cos(vetor_angulo(2) - vetor_angulo(3))+ parte_imaginaria(2,3)*sin(vetor_angulo(2) - vetor_angulo(3)))+vetor_tensao(2)*vetor_tensao(4)*(parte_real(2,4)*cos(vetor_angulo(2) - vetor_angulo(4))+ parte_imaginaria(2,4)*sin(vetor_angulo(2) - vetor_angulo(4)));
p3_calc =  vetor_tensao(3)*vetor_tensao(1)*(parte_real(3,1)*cos(vetor_angulo(3) - vetor_angulo(1))+ parte_imaginaria(3,1)*sin(vetor_angulo(3) - vetor_angulo(1)))+vetor_tensao(3)*vetor_tensao(2)*(parte_real(3,2)*cos(vetor_angulo(3) - vetor_angulo(2))+ parte_imaginaria(3,2)*sin(vetor_angulo(3) - vetor_angulo(2)))+vetor_tensao(3)*vetor_tensao(3)*(parte_real(3,3)*cos(vetor_angulo(3) - vetor_angulo(3))+ parte_imaginaria(3,3)*sin(vetor_angulo(3) - vetor_angulo(3)))+vetor_tensao(3)*vetor_tensao(4)*(parte_real(3,4)*cos(vetor_angulo(3) - vetor_angulo(4))+ parte_imaginaria(3,4)*sin(vetor_angulo(3) - vetor_angulo(4)));
p4_calc =  vetor_tensao(4)*vetor_tensao(1)*(parte_real(4,1)*cos(vetor_angulo(4) - vetor_angulo(1))+ parte_imaginaria(4,1)*sin(vetor_angulo(4) - vetor_angulo(1)))+vetor_tensao(4)*vetor_tensao(2)*(parte_real(4,2)*cos(vetor_angulo(4) - vetor_angulo(2))+ parte_imaginaria(4,2)*sin(vetor_angulo(4) - vetor_angulo(2)))+vetor_tensao(4)*vetor_tensao(3)*(parte_real(4,3)*cos(vetor_angulo(4) - vetor_angulo(3))+ parte_imaginaria(4,3)*sin(vetor_angulo(4) - vetor_angulo(3)))+vetor_tensao(1)*vetor_tensao(4)*(parte_real(4,4)*cos(vetor_angulo(4) - vetor_angulo(4))+ parte_imaginaria(4,4)*sin(vetor_angulo(4) - vetor_angulo(4)));
%Potencia REATIVA
q1_calc =  vetor_tensao(1)*vetor_tensao(1)*(parte_imaginaria(1,1)*cos(vetor_angulo(1) - vetor_angulo(1))- parte_real(1,1)*sin(vetor_angulo(1) - vetor_angulo(1)))+vetor_tensao(1)*vetor_tensao(2)*(parte_imaginaria(1,2)*cos(vetor_angulo(1) - vetor_angulo(2))- parte_real(1,2)*sin(vetor_angulo(1) - vetor_angulo(2)))+vetor_tensao(1)*vetor_tensao(3)*(parte_imaginaria(1,3)*cos(vetor_angulo(1) - vetor_angulo(3))- parte_real(1,3)*sin(vetor_angulo(1) - vetor_angulo(3)))+vetor_tensao(1)*vetor_tensao(4)*(parte_imaginaria(1,4)*cos(vetor_angulo(1) - vetor_angulo(4))- parte_real(1,4)*sin(vetor_angulo(1) - vetor_angulo(4)));
q2_calc =  vetor_tensao(2)*vetor_tensao(1)*(parte_imaginaria(2,1)*cos(vetor_angulo(2) - vetor_angulo(1))- parte_real(2,1)*sin(vetor_angulo(2) - vetor_angulo(1)))+vetor_tensao(2)*vetor_tensao(2)*(parte_imaginaria(2,2)*cos(vetor_angulo(2) - vetor_angulo(2))- parte_real(2,2)*sin(vetor_angulo(2) - vetor_angulo(2)))+vetor_tensao(2)*vetor_tensao(3)*(parte_imaginaria(2,3)*cos(vetor_angulo(2) - vetor_angulo(3))- parte_real(2,3)*sin(vetor_angulo(2) - vetor_angulo(3)))+vetor_tensao(2)*vetor_tensao(4)*(parte_imaginaria(2,4)*cos(vetor_angulo(2) - vetor_angulo(4))- parte_real(2,4)*sin(vetor_angulo(2) - vetor_angulo(4)));
q3_calc =  vetor_tensao(3)*vetor_tensao(1)*(parte_imaginaria(3,1)*cos(vetor_angulo(3) - vetor_angulo(1))- parte_real(3,1)*sin(vetor_angulo(3) - vetor_angulo(1)))+vetor_tensao(3)*vetor_tensao(2)*(parte_imaginaria(3,2)*cos(vetor_angulo(3) - vetor_angulo(2))- parte_real(3,2)*sin(vetor_angulo(3) - vetor_angulo(2)))+vetor_tensao(3)*vetor_tensao(3)*(parte_imaginaria(3,3)*cos(vetor_angulo(3) - vetor_angulo(3))- parte_real(3,3)*sin(vetor_angulo(3) - vetor_angulo(3)))+vetor_tensao(3)*vetor_tensao(4)*(parte_imaginaria(3,4)*cos(vetor_angulo(3) - vetor_angulo(4))- parte_real(3,4)*sin(vetor_angulo(3) - vetor_angulo(4)));
q4_calc =  vetor_tensao(4)*vetor_tensao(1)*(parte_imaginaria(4,1)*cos(vetor_angulo(4) - vetor_angulo(1))- parte_real(4,1)*sin(vetor_angulo(4) - vetor_angulo(1)))+vetor_tensao(4)*vetor_tensao(2)*(parte_imaginaria(4,2)*cos(vetor_angulo(4) - vetor_angulo(2))- parte_real(4,2)*sin(vetor_angulo(4) - vetor_angulo(2)))+vetor_tensao(4)*vetor_tensao(3)*(parte_imaginaria(4,3)*cos(vetor_angulo(4) - vetor_angulo(3))- parte_real(4,3)*sin(vetor_angulo(4) - vetor_angulo(3)))+vetor_tensao(1)*vetor_tensao(4)*(parte_imaginaria(4,4)*cos(vetor_angulo(4) - vetor_angulo(4))- parte_real(4,4)*sin(vetor_angulo(4) - vetor_angulo(4)));

