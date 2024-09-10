%Método de resolução. O script precisa falar quais são os valores de A e b.
%Como exemplo fiz A = A1 e b = b1
A = A2;
b = b2;
mat_conc = [A,b];
mat_conc;
[linha,coluna] = size(mat_conc);
for k =1: (linha-1);
    for i = (k+1): linha;
        pivo = mat_conc(i,k)/mat_conc(k,k);
        for j = (k):(linha+1);
            mat_conc(i,j) = mat_conc(i,j) - pivo*mat_conc(k,j);
            mat_conc;
        end
    end
end
mat_conc;
solucao = mat_conc(:,end);
for k = linha:-1 : 1
    solucao(k) = mat_conc(k,end)
    solucao(k)
    for i = (k+1):linha
        solucao(k) = solucao(k)- mat_conc(k,i)*solucao(i)
    end
    solucao(k) = solucao(k)/mat_conc(k,k)
end
%Meu vetor solução é solucao abaixo
solucao


