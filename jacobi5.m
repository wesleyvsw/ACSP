%Método de resolução. O script precisa falar quais são os valores de A e b.
%Como exemplo fiz A = A1 e b = b1
A  = A2;
b = b2;
[linha,coluna] = size(A);
%O algoritmo divide  b- Ax por cada elemento da diagonal principal. Posso
%resolver o prblema de forma matricial
%Criando uma matriz com os elementos da diagonal principal
diagonal = A.*eye(linha);
d2 =diagonal;
%Invertendo essa matriz diagonal sem usar a função inv
for k= 1:linha
    d2(k,k) = 1/d2(k,k);
end
%Criando a matriz B que é uma matriz com diagonal igual a zero, com todos
%os termos dividido por cada elemento de sua respectiva diagonal principal
B= eye(linha) - d2*A;
%Dividindo o vetor b por cada elemento elemento diagonal equivalente da
%matriz A
g = d2*b;
x1 = zeros(linha,1);
contador =1;
erro =10000;
while erro > 0.01 & contador <100 
    %equacao de atualizacao
    x2 = B*x1 +g;
    erro = max(abs(x2-x1))/max(abs(x2));
    contador = contador +1;
    x1 = x2;
end

%Meu vetor solução é x2
x2

