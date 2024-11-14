clc
clear
l = 2 ;
m = 100 ;
dx = l/m ; 
alpha = 11.234*power(10,-5) ; %pure_copper
k = 386 ; %pure_copper
dt = 0.6 ; 
beta = (alpha*dt)/(power(dx,2)) ; 
Q = 50000;%ConstantFlux
Tw = 400 ; 
To = 800 ; 
n_max = 800 ; 
T = zeros(m+1 , n_max) ; 
T(: , 1) = To ;
T(1 , :) = Tw ; 
Q_transfer = zeros(1 , n_max) ; 
B = zeros(m-1) ;
TK = zeros(m+1) ;
TK_1 = zeros(m+1) ; 
TK = T ; 
TK_1 = T ;
eps = 0.001 ; 

for n=1 : n_max
    for i=2 : m-1
        B(i-1) = -beta*T(i+1 , n) + (2*beta-1)*T(i , n) - beta*T(i-1 , n) ;
    end
    B(m-1) = -beta*(dx*Q/k) + (beta-1)*T(m,n) - beta*T(m-1,n) ; 

    error = 10 ;
    num_loop = 0 ;

    while error > eps 

        num_loop = num_loop + 1 ;

        for i=2 : m-1
            TK_1(i) = 1/(2*beta+1)*(-B(i-1) + beta*TK_1(i+1) + beta*TK_1(i-1)) ;
        end

        TK_1(m) = (1/(beta+1))*(-B(m-1) + beta*(dx*Q/k) + beta*TK_1(m-1)) ;
        TK_1(m+1) = (dx*Q/k) + TK_1(m) ;

        error = norm(TK_1 - TK) ; 
        TK = TK_1 ;

        for i=1 : m+1
            T(i , n+1) = TK_1(i) ; 
        end

        num_loop
        n
        

    end    

         Q_transfer(1,n) = k*(T(2,n) - T(1,n))/dx ; %Output flux from the left side of the rod

     
end 



for n=1 :80: n_max 
    hold on ;
    plot(T(: , n)) ; 
    title('Temperature changes over time');
end

figure ;
hold on ;
plot(Q_transfer) ; 
title('heat flux');
