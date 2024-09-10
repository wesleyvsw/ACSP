% Definindo a matriz A e o vetor b
A = A1;  
b = b1;  

contador = 0;         
max_cont = 5000;     
erro = 0.0001;       

[linha, coluna] = size(A);  
C = A;                      

x1_gs = zeros(coluna, 1);   
for i = 1:coluna
    C(i, i) = 0;            
end

d = zeros(coluna, 1);       
for i = 1:coluna
    C(i, 1:coluna) = C(i, 1:coluna) / A(i, i);  
    d(i) = b(i) / A(i, i);                      
end

while (1)
    x2_gs = x1_gs;  
    for i = 1:coluna
        x1_gs(i) = d(i) - C(i, :) * x1_gs;  
        if x1_gs(i) ~= 0
            ea(i) = abs((x1_gs(i) - x2_gs(i)) / x1_gs(i));  
        else
            ea(i) = 0;  
        end
    end
    contador = contador + 1; 
    if max(ea) <= erro || contador >= max_cont 
        break;
    end
end

% Solução obtida pelo Método de Gauss-Seidel
x1_gs