A  = [15 5 -5;4 10 1; 2 -2 8];
b = [30 23 -10]';
loop = 0;
max_loop = 100;
erro = 0.001;
[linha,coluna] = size(A);
mat_aux = A
%Zerando elementos da diagonal principal e criando o vetor solução com
%valores iniciais iguais a zero
for k = 1:linha
    mat_aux(k,k) = 0;
    x(i) = 0
end
%invertendo o vetor solução
x = x'
for k = 1:linha
    mat_aux(k,1:linha) = mat_aux(k,1:linha)/A(k,k);
    d(i) = b(k)/A(k,k);
end
while (1)
    x0 = x;
    for i =1:linha
        x(i) = d(i)- mat_aux(i,:)*x;
        if x(i) ~= 0
            erro_calculado(i) = abs((x(i)-x0(i))/x(i));
        end
    end
    loop = loop + 1
    if max(erro_calculado)<=erro ||loop>=max_loop,break,end
end




