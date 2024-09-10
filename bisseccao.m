x0 = 0
x1 = 1
funcao = [2,3,5]
raiz = 'aaa'
iter = 0
erro = 0.0001
while polyval(ponto_medio)>erro & iter<100
    ponto_medio = (x0+x1)/2
    if polyval(ponto_medio) <= erro
        raiz = ponto_medio
        else
        if polyval(funcao,ponto_medio)<0
            x0 = valor_medio
        else
            x1 = valor_medio
        end
    end
    iter = iter+1
end