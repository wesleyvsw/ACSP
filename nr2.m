matriz_impedancias = [0   0+0.015i  0   0;0.015i 0 0.0009+0.014i 0; 0 0.0009+0.014i 0 0.015i; 0 0 0.015i 0];
matriz_adm_imp = 1./matriz_impedancias;%matriz com o inverso de cada impedancia
num_barras = 4;
contador = 0
b_shunt = 1.2i;
potencias_especificadas = [4;0;-2;0;-1]%[p1,p2,p3,q2,q3]
vetor_tensao = [1;1;1;1];
vetor_angulo = [0;0;0;0];
vetor_v_theta = [vetor_angulo(1);vetor_angulo(2);vetor_angulo(3);vetor_tensao(2);vetor_tensao(3)];%[theta1,theta2,theta3,v2,v3]
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
p_calc = zeros(num_barras, 1)
'for i = 1:num_barras
    % Cálculo da potência ativa para a barra i
    for j = 1:num_barras
      p_calc(i) = p_calc(i)+ vetor_tensao(i)vetor_tensao(j) * (parte_real(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)) + parte_imaginaria(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)));
    end
end'
'%Potencias reativas
q_calc=zeros(num_barras,1);
for i = 1:num_barras
    % Cálculo da potência ativa para a barra i
    for j = 1:num_barras
      q_calc(i) = q_calc(i)+ vetor_tensao(i)vetor_tensao(j) * (parte_real(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)) - parte_imaginaria(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)));
    end
end'
'potencias_calculadas =[p_calc(1);p_calc(2);p_calc(3);q_calc(2);p_calc(3)];% Vetor de potencias calculadas para realizar a atualização
########################################################################################
#####While principal do processo de Newton Rapshon######################################
########################################################################################
while max(abs(potencias_especificadas-potencias_calculadas))>0.001 && contador < 100
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Matriz Jacobiana
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Matriz H - dP/dTheta
  j11 = -1*vetor_tensao(1)*vetor_tensao(1)*parte_imaginaria(1,1)-q_calc(1);
  j22 = -1*vetor_tensao(2)*vetor_tensao(2)*parte_imaginaria(2,2)-q_calc(2);
  j33 = -1*vetor_tensao(3)*vetor_tensao(3)*parte_imaginaria(3,3)-q_calc(3);
  j12 = vetor_tensao(1)*vetor_tensao(2)*(parte_real(1,2)*sin(vetor_angulo(1)-vetor_angulo(2))-parte_imaginaria(1,2)*cos(vetor_angulo(1)-vetor_angulo(2)));
  j13 = vetor_tensao(1)*vetor_tensao(3)*(parte_real(1,3)*sin(vetor_angulo(1)-vetor_angulo(3))-parte_imaginaria(1,3)*cos(vetor_angulo(1)-vetor_angulo(3)));
  j21 = vetor_tensao(2)*vetor_tensao(1)*(parte_real(2,1)*sin(vetor_angulo(2)-vetor_angulo(1))-parte_imaginaria(2,1)*cos(vetor_angulo(2)-vetor_angulo(1)));
  j23 = vetor_tensao(2)*vetor_tensao(3)*(parte_real(2,3)*sin(vetor_angulo(2)-vetor_angulo(3))-parte_imaginaria(2,3)*cos(vetor_angulo(2)-vetor_angulo(3)));
  j31 = vetor_tensao(3)*vetor_tensao(1)*(parte_real(3,1)*sin(vetor_angulo(3)-vetor_angulo(1))-parte_imaginaria(3,1)*cos(vetor_angulo(3)-vetor_angulo(1)));
  j32 = vetor_tensao(3)*vetor_tensao(2)*(parte_real(3,2)*sin(vetor_angulo(3)-vetor_angulo(2))-parte_imaginaria(3,2)*cos(vetor_angulo(3)-vetor_angulo(2)));
  %Matriz N - dP/dV
  j24 = (p_calc(2)+vetor_tensao(2)*vetor_tensao(2)*parte_real(2,2))/vetor_tensao(2);
  j35 = (p_calc(3)+vetor_tensao(3)*vetor_tensao(3)*parte_real(3,3))/vetor_tensao(3);
  j14 = vetor_tensao(1)*(parte_real(1,2)*cos(vetor_angulo(1)-vetor_angulo(2))+parte_imaginaria(1,2)*sin(vetor_angulo(1)-vetor_angulo(2)));
  j15 = vetor_tensao(1)*(parte_real(1,3)*cos(vetor_angulo(1)-vetor_angulo(3))+parte_imaginaria(1,3)*sin(vetor_angulo(1)-vetor_angulo(3)));
  j25 = vetor_tensao(2)*(parte_real(2,3)*cos(vetor_angulo(2)-vetor_angulo(3))+parte_imaginaria(2,3)*sin(vetor_angulo(2)-vetor_angulo(3)));
  j34 = vetor_tensao(3)*(parte_real(3,2)*cos(vetor_angulo(3)-vetor_angulo(2))+parte_imaginaria(3,2)*sin(vetor_angulo(3)-vetor_angulo(2)));
  %MATRIZ M
  j42 = p_calc(2)- vetor_tensao(2)*vetor_tensao(2)*parte_real(2,2);
  j53 = p_calc(3)- vetor_tensao(3)*vetor_tensao(3)*parte_real(3,3);
  j41 = -1*vetor_tensao(2)*vetor_tensao(1)*(parte_real(2,1)*cos(vetor_angulo(2)-vetor_angulo(1))+parte_imaginaria(2,1)*sin(vetor_angulo(2)-vetor_angulo(1)));
  j43 = -1*vetor_tensao(2)*vetor_tensao(3)*(parte_real(2,3)*cos(vetor_angulo(2)-vetor_angulo(3))+parte_imaginaria(2,3)*sin(vetor_angulo(2)-vetor_angulo(3)));
  j51 = -1*vetor_tensao(3)*vetor_tensao(1)*(parte_real(3,1)*cos(vetor_angulo(3)-vetor_angulo(1))+parte_imaginaria(3,1)*sin(vetor_angulo(3)-vetor_angulo(1)));
  j52 = -1*vetor_tensao(3)*vetor_tensao(2)*(parte_real(3,2)*cos(vetor_angulo(3)-vetor_angulo(2))+parte_imaginaria(3,2)*sin(vetor_angulo(3)-vetor_angulo(2)));
  %MATRIZ L
  j44 = (q_calc(2)-vetor_tensao(2)*vetor_tensao(2)*parte_imaginaria(2,2))/vetor_tensao(2);
  j55 = (q_calc(3)-vetor_tensao(3)*vetor_tensao(3)*parte_imaginaria(3,3))/vetor_tensao(3);
  j45 = vetor_tensao(2)*(parte_real(2,3)*sin(vetor_angulo(2)-vetor_angulo(3))-parte_imaginaria(2,3)*cos(vetor_angulo(2)-vetor_angulo(3)));
  j54 = vetor_tensao(3)*(parte_real(3,2)*sin(vetor_angulo(3)-vetor_angulo(2))-parte_imaginaria(3,2)*cos(vetor_angulo(3)-vetor_angulo(2)));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  matriz_jacobiana = [j11,j12,j13,j14,j15;j21,j22,j23,j24,j25;j31,j32,j33,j34,j35;j41,j42,j43,j44,j45;j51,j52,j53,j54,j55];
  vetor_v_theta = vetor_v_theta + inv(matriz_jacobiana)*(potencias_especificadas-potencias_calculadas);%atualizando o vetor theta v
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%Atualizando o vetor tensao e angulo com os elementos do vetor theta_v%%%%%%%%%%%%%%%%%%%%%
  vetor_tensao(2) = vetor_v_theta(4);
  vetor_tensao(3) = vetor_v_theta(5);
  vetor_angulo(1) = vetor_v_theta(1);
  vetor_angulo(2) = vetor_v_theta(2);
  vetor_angulo(3) = vetor_v_theta(3);
  for i = 1:num_barras
    % Cálculo da potência ativa para a barra i
    for j = 1:num_barras
      p_calc(i) = p_calc(i)+ vetor_tensao(i)vetor_tensao(j) * (parte_real(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)) + parte_imaginaria(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)));
    end
  end
  %Potencias reativas
  for i = 1:num_barras
    % Cálculo da potência ativa para a barra i
    for j = 1:num_barras
      q_calc(i) = q_calc(i)+ vetor_tensao(i)vetor_tensao(j) * (parte_real(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)) - parte_imaginaria(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)));
    end
  end
  potencias_calculadas =[p_calc(1);p_calc(2);p_calc(3);q_calc(2);p_calc(3)];% Vetor de potencias calculadas para realizar a atualização
  contador = contador+1
 end'
