clear all
close all

b = 5; % in m
m = 1;
Q = 35; % in m^3/sec
n = 0.013;
So = 0.004;
g = 9.81;
y_given = 3.4;

yn = trapz_normal_depth(Q,b,m,n,So);
yc = trapz_critical_depth(Q,b,m);

y = y_given;
A = y*b + y^2*m ;
P = b + 2*(y^2+(y*m)^2)^0.5;
Rh = A/P;
V = Q/A;
vel_head = V^2/2/g;
E = vel_head+y;
Se = n^2*V^2/Rh^(3/4);
L = 0;

wanted_diff = 1000;
del_y = (y-yc)/wanted_diff;

for i = y-del_y:-del_y:yc
    y = [y, i];
    A = [A, i*b+i^2*m] ;
    P = [P, b+2*(i^2+(i*m)^2)^0.5] ;
    Rh = [Rh, A(length(A))/P(length(P))];
    V = [V, Q/A(length(A))];
    vel_head = [vel_head, V(length(V))^2/2/g];
    E = [E, vel_head(length(vel_head)) + i ];
    Se = [Se, n^2*V(length(V))^2/Rh(length(Rh))^(3/4)];
    del_L = (E(length(E)-1)-E(length(E)))/(So-(abs(Se(length(Se))-Se(length(Se)-1)))/2);
    L = [L, L(length(L))+del_L];
end

fprintf('Distance between dam and the hydraulic jump L = %.5f \n', L(length(L)));

function result = trapz_critical_depth(Q,b,m);
    g = 9.81;    
    y = [0,100];    
    tolerance = 1;    
    rhs = Q^2/g;

    while tolerance>10^-5
        ym = (y(1)+y(2))/2;
        T = 2*ym*m + b;
        A = b*ym + ym^2*m;
        lhs = A^3/T;        
        tolerance = abs(y(1)-y(2));
    
        if lhs>rhs
            y =[y(1),ym];
        else
            y =[ym,y(2)];
        end
    end
    result = (y(1)+y(2))/2;
end



function result = trapz_normal_depth(Q,b,m,n,So);
    y = [0,100];    
    tolerance = 1;
    rhs = Q*n/So^0.5;

    while tolerance > 10^-5    
        
        ym = (y(1)+y(2))/2;    
        A = ym*b + ym^2*m ;
        P = b + 2*(ym^2+(ym*m)^2)^0.5 ;        
        lhs = A*(A/P)^(2/3);        
        tolerance = abs(y(1)-y(2));
    
        if lhs>rhs
            y =[y(1),ym];
        else
            y =[ym,y(2)];
        end
    end
    result = (y(1)+y(2))/2;
end