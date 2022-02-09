close all
clear all

Za = 100; %in m
Zb = 110; %in m

e = 0.26/1000; %in m
D = 0.13; %in m
L = 20; %in m
vis = 10^-6; %in m^2/s
g = 9.81;
%pump characteristic curve

del_press = [3.90 3.20 2.55 2.00 1.20]; %in bar
in_press =  [-1.20 -1.80 -2.10 -2.30 -2.40]; %in bar

Q_in_l = [50 100 150 200 250]; % in l/s
Q = Q_in_l*0.001; % in m^3/s

% 1 bar = 10.1974 metres of water in 20°C

pump_press_in_bar = del_press-in_press; % in bar
pump_press = pump_press_in_bar.*10.1974; % in m

plot(Q,pump_press);

%system curve

area = pi*D^2/4; % in m^2

vel = Q./area; % in m^3/s/m^2 = m/s
RN = vel.*D/vis;

%f(i)=0.25./((log10(ed./3.7+5.74./(Re(i).^0.9))).^2);
f = zeros(1,length(Q));

for i = 1:1:length(Q)
    f(i) = 0.25/( (log10( (e/D/3.7) + (5.74/RN(i)^0.9) ) )^2 );
end

for i = 1:1:length(Q);
    hl(i) = f(i)*(L/D)*(vel(i)^2/2/g);
end

for i = 1:1:length(Q);
    hp(i) = (Zb-Za)+hl(i);
end

hold on
plot(Q,hp);xlabel("Q (m^3/s)");ylabel("Head (m)");


%intersection
% we know that it is between Q = 0.15 and Q = 0.2 in pump characteristic
% curve and it is a line

m = (pump_press(3)-pump_press(4))/(Q(3)-Q(4));
b = pump_press(3)-m*Q(3);

%eqn is y = mx+b

Q_pred =[0.15 0.2];
tolerance = 100;
while tolerance > 10^-5
    
    Q_pred_mid = abs(Q_pred(1)+Q_pred(2))/2;
    
    hp_target = m*Q_pred_mid+b;

    vel_pred = Q_pred_mid/area;
    Rn_pred = vel_pred*D/vis;
    f_pred = 0.25/( (log10( (e/D/3.7) + (5.74/Rn_pred^0.9) ) )^2 );
    hl_pred = f_pred*L/D*vel_pred^2/2/g;
    hp_pred = (Zb-Za)+hl_pred;

    tolerance = abs(hp_target-hp_pred);
    if hp_target>hp_pred
        Q_pred =[Q_pred_mid Q_pred(2)];
    else
        Q_pred =[Q_pred(1) Q_pred_mid];
    end
end

fprintf("The system Flow Rate in m^3/s is : %.4f and the required head in meter is : %.3f \n",[Q_pred_mid,hp_pred]);


