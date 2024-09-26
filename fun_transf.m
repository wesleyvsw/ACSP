%Função step no octave
k = 1.2;
tau = 0.1;
tmax = 10.0;
num = k;
den = [tau 1];
fun_transf = tf(num,den)
vetor_tempo = [0:0.001:tmax];
[y,y1] = step(fun_transf,vetor_tempo);
plot(vetor_tempo,y)

