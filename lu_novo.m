A  = [4 -9 2;2 -4 4; -1 2 2];
b = [2 3 1]';
mat_conc = [A,b];
mat_conc;
[linha,coluna] = size(mat_conc);
mat_l = eye(linha)
for k =1: (linha-1);
    for i = (k+1): linha;
        pivo = mat_conc(i,k)/mat_conc(k,k);
        mat_l(i,k) = pivo;
        for j = (k):(linha+1);
            mat_conc(i,j) = mat_conc(i,j) - pivo*mat_conc(k,j);
            mat_conc;
        end
    end
end
mat_conc;
mat_conc = mat_conc(:, 1:end-1)
mat_l
solucao_y = b;
for k = 2 : linha
    solucao_y(k) = solucao_y(k)-(mat_l(k,:)*solucao_y-solucao_y(k));
    solucao_y(k);
end
mat_conc = [mat_conc,solucao_y]
solucao = mat_conc(:,end);
for k = linha:-1 : 1
    solucao(k) = mat_conc(k,end)
    solucao(k)
    for i = (k+1):linha
        solucao(k) = solucao(k)- mat_conc(k,i)*solucao(i)
    end
    solucao(k) = solucao(k)/mat_conc(k,k)
end
