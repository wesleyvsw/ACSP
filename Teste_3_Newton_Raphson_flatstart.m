% Definir a matriz de impedâncias
matriz_impedancias = [0   0+0.015i  0   0; 0+0.015i 0 0.0009+0.014i 0; 0 0.0009+0.014i 0 0+0.015i; 0 0 0+0.015i 0];

% Calcular a matriz de admitâncias invertendo cada impedância
matriz_adm_imp = 1./matriz_impedancias;

num_barras = 4;
contador = 0;
b_shunt = 1.2i;

% Potências especificadas [p1, p2, p3, q2, q3]
potencias_especificadas = [4; 0; -2; 0; -1];

% Tensões iniciais e ângulos iniciais
vetor_tensao = [1; 1; 1; 1];
vetor_angulo = [0; 0; 0; 0];
vetor_v_theta = [vetor_angulo(1); vetor_angulo(2); vetor_angulo(3); vetor_tensao(2); vetor_tensao(3)];

% Tratar infinitos e NaNs na matriz de admitâncias
mascara_infinitos = isinf(matriz_adm_imp);
matriz_adm_imp = matriz_adm_imp .* ~mascara_infinitos;
mascara_nan = isnan(matriz_adm_imp);
matriz_adm_imp(mascara_nan) = 0;

% Criar a matriz de admitâncias
mat_adm = matriz_adm_imp;
mat_adm(2,3) = 2 * mat_adm(2,3);
mat_adm(3,2) = 2 * mat_adm(3,2); % Linhas em paralelo

% Adicionar elementos shunt nas barras 2 e 3
mat_adm(2,2) = mat_adm(2,2) + b_shunt;
mat_adm(3,3) = mat_adm(3,3) + b_shunt;

% Atualizar os elementos diagonais
[linhas, colunas] = size(mat_adm);
for i = 1:linhas
    soma_linha = sum(mat_adm(i, :));
    mat_adm(i, i) = soma_linha;
end

% Multiplicar elementos fora da diagonal principal por -1
for x = 1:linhas
    for y = 1:colunas
        if x ~= y
            mat_adm(x, y) = mat_adm(x, y) * -1;
        end
    end
end

% Partes real e imaginária da matriz de admitâncias
parte_real = real(mat_adm);
parte_imaginaria = imag(mat_adm);

%%%%%%%%%%%%%%%%%%%%%%%% MÉTODO DE NEWTON %%%%%%%%%%%%%%%%%%%%%%

% Calcular potências ativas iniciais
p_calc = zeros(num_barras, 1);
for i = 1:num_barras
    for j = 1:num_barras
        p_calc(i) = p_calc(i) + vetor_tensao(i) * vetor_tensao(j) * (parte_real(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)) + parte_imaginaria(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)));
    end
end

% Calcular potências reativas iniciais
q_calc = zeros(num_barras, 1);
for i = 1:num_barras
    for j = 1:num_barras
        q_calc(i) = q_calc(i) + vetor_tensao(i) * vetor_tensao(j) * (parte_real(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)) - parte_imaginaria(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)));
    end
end

% Potências calculadas iniciais para a atualização no método de Newton-Raphson
potencias_calculadas = [p_calc(1); p_calc(2); p_calc(3); q_calc(2); q_calc(3)];

########################################################################################
##### Loop principal do processo de Newton-Raphson #####################################
########################################################################################

while max(abs(potencias_especificadas - potencias_calculadas)) > 0.001 && contador < 100
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matriz Jacobiana
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matriz H - dP/dTheta
    j11 = -vetor_tensao(1) * vetor_tensao(1) * parte_imaginaria(1,1) - q_calc(1);
    j22 = -vetor_tensao(2) * vetor_tensao(2) * parte_imaginaria(2,2) - q_calc(2);
    j33 = -vetor_tensao(3) * vetor_tensao(3) * parte_imaginaria(3,3) - q_calc(3);
    j12 = vetor_tensao(1) * vetor_tensao(2) * (parte_real(1,2) * sin(vetor_angulo(1) - vetor_angulo(2)) - parte_imaginaria(1,2) * cos(vetor_angulo(1) - vetor_angulo(2)));
    j13 = vetor_tensao(1) * vetor_tensao(3) * (parte_real(1,3) * sin(vetor_angulo(1) - vetor_angulo(3)) - parte_imaginaria(1,3) * cos(vetor_angulo(1) - vetor_angulo(3)));
    j21 = vetor_tensao(2) * vetor_tensao(1) * (parte_real(2,1) * sin(vetor_angulo(2) - vetor_angulo(1)) - parte_imaginaria(2,1) * cos(vetor_angulo(2) - vetor_angulo(1)));
    j23 = vetor_tensao(2) * vetor_tensao(3) * (parte_real(2,3) * sin(vetor_angulo(2) - vetor_angulo(3)) - parte_imaginaria(2,3) * cos(vetor_angulo(2) - vetor_angulo(3)));
    j31 = vetor_tensao(3) * vetor_tensao(1) * (parte_real(3,1) * sin(vetor_angulo(3) - vetor_angulo(1)) - parte_imaginaria(3,1) * cos(vetor_angulo(3) - vetor_angulo(1)));
    j32 = vetor_tensao(3) * vetor_tensao(2) * (parte_real(3,2) * sin(vetor_angulo(3) - vetor_angulo(2)) - parte_imaginaria(3,2) * cos(vetor_angulo(3) - vetor_angulo(2)));

    % Matriz N - dP/dV
    j14 = vetor_tensao(1) * (parte_real(1,2) * cos(vetor_angulo(1) - vetor_angulo(2)) + parte_imaginaria(1,2) * sin(vetor_angulo(1) - vetor_angulo(2)));
    j15 = vetor_tensao(1) * (parte_real(1,3) * cos(vetor_angulo(1) - vetor_angulo(3)) + parte_imaginaria(1,3) * sin(vetor_angulo(1) - vetor_angulo(3)));
    j24 = (p_calc(2) + vetor_tensao(2) * vetor_tensao(2) * parte_real(2,2)) / vetor_tensao(2);
    j25 = vetor_tensao(2) * (parte_real(2,3) * cos(vetor_angulo(2) - vetor_angulo(3)) + parte_imaginaria(2,3) * sin(vetor_angulo(2) - vetor_angulo(3)));
    j34 = vetor_tensao(3) * (parte_real(3,2) * cos(vetor_angulo(3) - vetor_angulo(2)) + parte_imaginaria(3,2) * sin(vetor_angulo(3) - vetor_angulo(2)));
    j35 = (p_calc(3) + vetor_tensao(3) * vetor_tensao(3) * parte_real(3,3)) / vetor_tensao(3);

    % Matriz M - dQ/dTheta
    j41 = -vetor_tensao(2) * vetor_tensao(1) * (parte_real(2,1) * cos(vetor_angulo(2) - vetor_angulo(1)) + parte_imaginaria(2,1) * sin(vetor_angulo(2) - vetor_angulo(1)));
    j42 = -vetor_tensao(2) * vetor_tensao(2) * parte_real(2,2) + p_calc(2);
    j43 = -vetor_tensao(2) * vetor_tensao(3) * (parte_real(2,3) * cos(vetor_angulo(2) - vetor_angulo(3)) + parte_imaginaria(2,3) * sin(vetor_angulo(2) - vetor_angulo(3)));
    j51 = -vetor_tensao(3) * vetor_tensao(1) * (parte_real(3,1) * cos(vetor_angulo(3) - vetor_angulo(1)) + parte_imaginaria(3,1) * sin(vetor_angulo(3) - vetor_angulo(1)));
    j52 = -vetor_tensao(3) * vetor_tensao(2) * (parte_real(3,2) * cos(vetor_angulo(3) - vetor_angulo(2)) + parte_imaginaria(3,2) * sin(vetor_angulo(3) - vetor_angulo(2)));
    j53 = -vetor_tensao(3) * vetor_tensao(3) * parte_real(3,3) + p_calc(3);

    % Matriz L - dQ/dV
    j44 = (q_calc(2) - vetor_tensao(2) * vetor_tensao(2) * parte_imaginaria(2,2)) / vetor_tensao(2);
    j45 = vetor_tensao(2) * (parte_real(2,3) * sin(vetor_angulo(2) - vetor_angulo(3)) - parte_imaginaria(2,3) * cos(vetor_angulo(2) - vetor_angulo(3)));
    j54 = vetor_tensao(3) * (parte_real(3,2) * sin(vetor_angulo(3) - vetor_angulo(2)) - parte_imaginaria(3,2) * cos(vetor_angulo(3) - vetor_angulo(2)));
    j55 = (q_calc(3) - vetor_tensao(3) * vetor_tensao(3) * parte_imaginaria(3,3)) / vetor_tensao(3);

    % Matriz Jacobiana completa
    matriz_jacobiana = [j11, j12, j13, j14, j15;
                        j21, j22, j23, j24, j25;
                        j31, j32, j33, j34, j35;
                        j41, j42, j43, j44, j45;
                        j51, j52, j53, j54, j55];

    % Atualizar ângulos e magnitudes de tensão
    delta = matriz_jacobiana \ (potencias_especificadas - potencias_calculadas);
    vetor_v_theta = vetor_v_theta + delta;

    % Atualizar magnitudes de tensão e ângulos
    vetor_angulo(1) = vetor_v_theta(1);
    vetor_angulo(2) = vetor_v_theta(2);
    vetor_angulo(3) = vetor_v_theta(3);
    vetor_tensao(2) = vetor_v_theta(4);
    vetor_tensao(3) = vetor_v_theta(5);

    % Recalcular potências ativas e reativas
    p_calc = zeros(num_barras, 1);
    for i = 1:num_barras
        for j = 1:num_barras
            p_calc(i) = p_calc(i) + vetor_tensao(i) * vetor_tensao(j) * (parte_real(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)) + parte_imaginaria(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)));
        end
    end

    q_calc = zeros(num_barras, 1);
    for i = 1:num_barras
        for j = 1:num_barras
            q_calc(i) = q_calc(i) + vetor_tensao(i) * vetor_tensao(j) * (parte_real(i,j) * sin(vetor_angulo(i) - vetor_angulo(j)) - parte_imaginaria(i,j) * cos(vetor_angulo(i) - vetor_angulo(j)));
        end
    end

    potencias_calculadas = [p_calc(1); p_calc(2); p_calc(3); q_calc(2); q_calc(3)];
    contador = contador + 1;
end

% Converter ângulos de radianos para graus
vetor_angulo_deg = vetor_angulo * (180 / pi);

% Exibir os resultados
disp('Magnitudes de tensão e ângulos em cada barra:')
for i = 1:num_barras
    fprintf('Barra %d: V = %.4f pu, Ângulo = %.4f graus\n', i, vetor_tensao(i), vetor_angulo_deg(i));
end

disp('Potências ativa e reativa em cada barra:')
for i = 1:num_barras
    fprintf('Barra %d: P = %.4f pu, Q = %.4f pu\n', i, p_calc(i), q_calc(i));
end

fprintf('Número de iterações: %d\n', contador);
