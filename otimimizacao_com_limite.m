% Definindo as funções de custo
%pkg load optim %Adicionando a biblioteca
C1 = @(g1) 10 * g1^4 + 5;
C2 = @(g2) 5 * g2^4 + 15;

Aeq = [1, 1];
beq = 100;

%Função objetivo e limite das variaveis
funcao_objetivo = @(x) C1(x(1)) + C2(x(2));
lb = [0; 0];
ub = [40,80];

options = optimset('Display', 'iter'); % Exibe informações durante a otimização
[x, fval] = fmincon(funcao_objetivo, [25; 75], [], [], Aeq, beq, lb, ub, [], options);

%Resultados
g1 = x(1);
g2 = x(2);
custo_total = fval;

disp(['Gerador 1 (g1): ', num2str(g1)]);
disp(['Gerador 2 (g2): ', num2str(g2)]);





