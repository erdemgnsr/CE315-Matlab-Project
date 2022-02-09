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

if So>0
    if yn>yc
        answer = "M";
    elseif yn==yc
        answer = "C";
    else
        answer = "S";
    end
elseif So==0
    answer ="H";
else 
    answer = "A";
end

if y_given>yn && y_given>yc
    answer = answer+"-1";
elseif (y_given>yn && y_given<yc) | (y_given<yn && y_given>yc)
    answer = answer+"-2";
else
    answer = answer+"-3";
end
fprintf("The flow classification result is:  ");
fprintf(answer+"\n");

dx = 50;
x = 0:dx:250;
len_x = length(x);

y = zeros(1,len_x); z = zeros(1,len_x); 
A = zeros(1,len_x); V = zeros(1,len_x); 
vel_head = zeros(1,len_x); P = zeros(1,len_x); 
Rh = zeros(1,len_x); Se_avg = zeros(1,len_x);
E = zeros(1,len_x) ;h_l = zeros(1,len_x - 1);

y(1) = y_given;
A(1) = y(1)*b + y(1)^2*m ;
z(1) = 0;
V(1) = Q/A(1);
vel_head(1) = V(1)^2/2/g;
P(1) = b + 2*(y(1)^2+(y(1)*m)^2)^0.5;
Rh(1) = A(1)/P(1);
Se(1) = n^2*V(1)^2/Rh(1)^(3/4);
E(1) = vel_head(1)+z(1)+y(1);


for i = 2:len_x
    tolerance = 10000;
    z(i) = x(i)*So;
    Y = [0,y(i-1)];
    while tolerance > 0.000001
        YM = (Y(2)+Y(1))/2;
        z(i) = x(i)*So;
        A(i) = YM*b + YM^2*m ;
        V(i) = Q/A(i);
        vel_head(i) = V(i)^2/2/g;
        P(i) = b + 2*(YM^2+(YM*m)^2)^0.5;
        Rh(i) = A(i)/P(i);
        Se(i) = n^2*V(i)^2/Rh(i)^(3/4);
        Se_avg = (Se(i)+Se(i-1))/2;
        delta_l = ( - (z(i-1) + vel_head(i-1)) + (z(i) + vel_head(i)) ) / So-Se_avg;
        h_l(i-1) = Se_avg*delta_l;
        E(i) = YM+z(i)+vel_head(i);
        tolerance = abs(E(i) - h_l(i-1) - E(i-1));
        if E(i) > h_l(i-1) + E(i-1)
            Y = [Y(1),YM];
        else
            Y = [YM,Y(2)];
        end
    end
    y(i) = YM;
    E(i-1) = E(i-1)+h_l(i-1);
end

for i = 1:length(y)
    fprintf("When x = %3d, the water surface profile y = %.3f\n",x(i),y(i));
end

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