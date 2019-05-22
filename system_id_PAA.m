clc;
clear all;
%% DATA
m1 = 300;
m2 = 60;
k1 = 16000;
k2 = 190000;
b1 = 1000;
%% Linearization into state space model
A = [0 1 0 0;-k1/m1 -b1/m1 k1/m1 b1/m1;0 0 0 1;k1/m2 b1/m2 (-k2-k1)/m2 -b1/m2];
B = [0;1/m1;0; -1/m2 ];
C = [1 0 0 0];
D = [0];
linsys_ss = ss(A,B,C,D);
%% converting state space into continuous transfer function
[num den] = ss2tf(A,B,C,D);
cont_tranf = tf(num,den)
%% continuous transfer function to discrete transfer function
Ts = .1;
disc_tranf = c2d(linsys_ss,Ts,'zoh');
Ad = disc_tranf.A;
Bd = disc_tranf.B;
Cd = disc_tranf.C;
Dd = disc_tranf.D;
[num1 den1] = ss2tf(Ad,Bd,Cd,Dd)
disc_z_inverse = filt([num1(1,1) num1(1,2) num1(1,3) num1(1,4)],[den1(1,1) den1(1,2) den1(1,3) den1(1,4) den1(1,5)])
%% Parametric adaptation algorithm with least square adaptation with forgetting factor
N = 6000;
t = 0:N-1;
%% Plant parameters
b1 = num1(1,1);
b2 = num1(1,2);
b3 = num1(1,3);
b4 = num1(1,4);
a1 = den1(1,2);
a2 = den1(1,3);
a3 = den1(1,4);
a4 = den1(1,5);
lamba = 0.99;
%% input and disturbance
w1 = 10;
w2 = 20;
w3 = 30;
u = sin(w1*t)+sin(w2*t)+sin(w3*t);
A = 0.05;
omega = A*randn(size(t));
%% Simulation of parametric Adaptation algorithm
f11 = 30;
f22 = 30;
f33 = 30;
f44 = 30;
f55 = 30;
f66 = 30; 
f77 = 30; 
f88 = 30;
F = diag([f11,f22,f33,f44,f55,f66,f77,f88]);

thetahat = zeros(8,N);
yv = zeros(1,N);
yhatv = zeros(1,N);
err = zeros(1,N);
f = zeros(1,N);

thhat = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';
y1=0;
y2=0;
y3=0;
y4=0;
u1=0;
u2=0;
u3=0;
u4=0;
omega1=0; 
phi=[0 0 0 0 0 0 0 0]';

for k=5:N 
    yv(k) = -a1*yv(k-1)-a2*yv(k-2)-a3*yv(k-3)-a4*yv(k-4)+b1*u(k-1)+b2*u(k-2)+b3*u(k-3)+b4*u(k-4)+omega(k);
    eo = yv(k)-thhat'*phi;
    F = (F-((F*(phi)*phi'*F)/(lamba+phi'*F*(phi))))*1/lamba;
    thhat = thhat+F*phi*eo;
    a1estimate(k) = thhat(1,:);
    a2estimate(k) = thhat(2,:);
    a3estimate(k) = thhat(3,:);
    a4estimate(k) = thhat(4,:);
    b1estimate(k) = thhat(5,:);
    b2estimate(k) = thhat(6,:);
    b3estimate(k) = thhat(7,:);
    b4estimate(k) = thhat(8,:);
    
    y4 = y3;
    y3 = y2;    
    y2 = y1;
    y1 = yv(k);
    u4 = u3;
    u3 = u2;
    u2 = u1;
    u1 = u(k);
    phi = [-y1 -y2 -y3 -y4 u1 u2 u3 u4]';
        % predictive error
    error(k) = eo;
end
figure(1)
plot(t,a1estimate,t,a1*ones(1,N))
title('a1 estimation')
figure(2)
plot(t,a2estimate,t,a2*ones(1,N))
title('a2 estimation')
figure(3)
plot(t,a3estimate,t,a3*ones(1,N))
title('a3 estimation')
figure(4)
plot(t,a4estimate,t,a4*ones(1,N))
title('a4 estimation')
figure(5)
plot(t,b1estimate,t,b1*ones(1,N))
title('b1 estimation')
figure(6)
plot(t,b2estimate,t,b2*ones(1,N))
title('b2 estimation')
figure(7)
plot(t,b3estimate,t,b3*ones(1,N))
title('b3 estimation')
figure(8)
plot(t,b4estimate,t,b4*ones(1,N))
title('b4 estimation')
figure(9)
plot(t,error,t,eo*ones(1,N))
title('error prediction')
