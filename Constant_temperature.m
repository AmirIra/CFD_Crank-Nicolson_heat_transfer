clc
clear
l = 2 ;
m = 100 ;
dx = l/m ; 
alpha = 11.234*power(10,-5) ; %pure_copper
k = 386 ; %pure_copper
dt = 0.5 ; 
beta = (alpha*dt)/(power(dx,2)) ; 
Tw = 400 ; 
Te = 400 ; 
To = 800 ; 
n_max = 600 ; 
T = zeros(m+1 , n_max) ; 
T(: , 1) = To ;
T(1 , :) = Tw ; 
T(m+1 , :) = Te ;
Q_transfer_w = zeros(1 , n_max) ; %Output flux from the left side of the rod
Q_transfer_e = zeros(1 , n_max) ; %The outlet flux is furnished from the right side
B = zeros(m-1) ;
TK = zeros(m+1) ;
TK_1 = zeros(m+1) ; 
TK = T ; 
TK_1 = T ;
eps = 0.0001 ; 

for n=1 : n_max
    for i=2 : m
        B(i-1) = -beta*T(i+1 , n) + (2*beta-1)*T(i , n) - beta*T(i-1 , n) ;
    end

    error = 10 ;
    num_loop = 0 ;

    while error > eps 

        num_loop = num_loop + 1 ;

        for i=2 : m
            TK_1(i) = 1/(2*beta+1)*(-B(i-1) + beta*TK_1(i+1) + beta*TK_1(i-1));
        end

        error = norm(TK_1 - TK) ; 
        TK = TK_1 ;

        for i=1 : m+1
            T(i , n+1) = TK_1(i) ; 
        end

        num_loop
        n
        

    end    


     Q_transfer_w(1,n) = abs(-k*(T(2,n)-Tw)/dx) ; %Output flux from the left side of the rod
     Q_transfer_e(1,n) = -abs(-k*(T(m,n)-Te)/dx) ; %Output flux from the right side of the rod 
end 

for n=1 :50: n_max 
    hold on ;
    plot(T(: , n)) ; 
    title('Temperature changes over time');
end
figure ;
hold on ;
plot(Q_transfer_e) ;
plot(Q_transfer_w) ;
% I have multiplied one of the graphs by the negative so that both graphs are on top of each other due to being the same
title('heat flux') ;




