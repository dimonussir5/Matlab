N = 25;
nodes = 30;
P = zeros(N+1, N+1); % создание матрицы из нулей
P(0+1,0+1) = 1; %P(0,0)
l0 = 9;
mu_pam = l0;
mu_per = 0.8; %интенсивность передачи пакетов по каналу
mu_pr = 15;% интенсивность обработки пакета в процессоре
L=9; %Число каналов по условию задачи
mu = zeros(1,nodes);
TO = 2/mu_per; %среднее время тайм оут
L_TO = 1/TO; % интенсивность

TU = 0.1 / mu_per; % время успешной доставки квитанции
L_TU = 1/TU;% интенсивность 
mu = [l0, mu_pr,mu_pam,mu_per,mu_per,mu_per,mu_per,mu_per,mu_per,mu_per,mu_per,mu_per,L_TU,L_TU,L_TU,L_TU,L_TU,L_TU,L_TU,L_TU,L_TU,L_TO,L_TO,L_TO,L_TO,L_TO,L_TO,L_TO,L_TO,L_TO, ];
m = [1,N,1,1,1,1,1,1,1,1,1,1,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N,N]; %число каналов
%расчет w
P_k = (1/L) * ones(1, L);
F = 0.1 * ones(1, L);
w_l = zeros(1, L);
w_m = zeros(1, L);
w_r = zeros(1, L);
w_k0 = 1;
w_l0 = 1;
w_m0 = 1;



for i = 1:L 
     w_l(i) = P_k(i) / (1-F(i));
     w_m(i) = P_k(i);
     w_r(i) = F(i) * P_k(i) / (1-F(i));
end

w = zeros(1,nodes);

w = [w_k0,w_l0,w_m0, w_l, w_m, w_r];
w_sum = sum(w);

P_0_N = zeros(1,nodes);
 %расчет по шагам
for i = 1:nodes
    mu_i = zeros(1,N);
    
    %расчет mu
    for n = 1:N
        if n < m(i)
             mu_i(n) = n * mu(i);
        else
            
             mu_i(n) = m(i) * mu(i); 
        end
    end
   
    t=zeros(1,N);
    l1 = zeros(1,N);
    for r = 1:N
        %шаг 1
       
    
        for n = 1:r
               t(r) = t(r) + n/mu_i(n) * P(n-1+1, r-1+1);
        end
    
        %шаг 2
        l1(r)= r / (w_sum * t(r) / w(1));
        %шаг 3 
        sum_P = 0;
        for n = 1:r
         %неизвестно какое w писать 
           P(n+1, r+1) = w(i) * l1(r) * P(n-1+1, r-1+1)/(w(1) * mu_i(n));
           sum_P = sum_P + P(n+1,r+1);
        end
        P(0+1,r+1) = 1-sum_P;
    
    end
    P_0_N(i) = P(0+1,N+1);
end
disp(P_0_N); % Вывод результата
n1 = N * (1-P_0_N(1));
l1 = l0*(1-P_0_N(1));
N_uk = N-n1;
t = N_uk / (l0 * (1-P_0_N(1)));
disp("a)вероятность отказа в приёме пакета в буферную память");
fprintf("P(0,N) =%.2f\n", P_0_N(1))
disp("б)среднее время пребывания пакета в буферной памяти");
fprintf("t = %.2f\n", t)
disp("в) среднее число пакетов в маршрутизаторе");
fprintf("N_uk = %.2f\n", N_uk)
disp("г) средняя интенсивность потока пакетов, занимающих буферную память маршрутизатора");
fprintf("l1 = %.2f\n", l1)


