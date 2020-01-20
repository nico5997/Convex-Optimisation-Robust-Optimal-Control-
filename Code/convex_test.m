clc;
close all
clear all

N = 35; % number of iterations
A1 = 2734e-6; %[m^2] Area of the cross-section of the first tank
A2 = 2734e-6; %[m^2] Area of the cross-section of the second tank
a1 = 7e-6; %[m^2] Area of of the cross-section of the first outlet
a2 = 7e-6; %[m^2] Area of of the cross-section of the second outlet
k_p = 50; %[V/m] Proportional constant between sensor and water height
k = 27e-7; %[m^3/Vs] Proportional contant of the pump
y_ss = 5; %[v] steady-state solution
u_ss = 1.82; %[V] steady-state solution
grv = 9.81; %[m.s^-2]
alpha1 = (-a1*k_p/A1)*sqrt(2*grv/k_p)*(1/2)*(1/sqrt(y_ss)); %[s^-1]
alpha12 = (a1*k_p/A2)*sqrt(2*grv/k_p)*(1/2)*(1/sqrt(y_ss)); %[s^-1]
alpha2 = (-a2*k_p/A2)*sqrt(2*grv/k_p)*(1/2)*(1/sqrt(y_ss)); %[s^-1]
beta = k*k_p/A1; %[S^-1]
A = [alpha1, 0; alpha12, alpha2]; %[S^-1]
B = [beta; 0]; %[S^-1]
d = [5; 5; 5; 5; 8.18; 1.82]; %[V]
C = [1 0 0;0 1 0;-1 0 0;0 -1 0; 0 0 0; 0 0 0];
D = [0 0 0 0 1 -1]';
Q = [0, 0, 0; 0, 2, 0; 0, 0, 0];
R = 0;
S = zeros(3,1);
T = 2; %[s] sampling time
L = 6; % Number of systems
sigma = 0.25; % relative standard deviation
u_0 = 4-1.82; %[V] initial condition
x_0 = [5-5; 1-5; 0]; %[V] initial condition

n = size(A);
n = n(1);
n = n+1;
m = size(B);
m = m(2);
p = length(d);

[F, phi_Q, g, C_p, d_p, x, gam] = OptRobustControl(A,B,C,D,d,Q,R,S,T,L,N,sigma,u_0,x_0);

%% Performance of robust control

t = (0:T:N*T); % Time vector

figure % Tanks level
for i = 1:3:L*n
hold on
plot(t,x(i:n*L+1:end)+y_ss);
end
for i = 2:3:L*n
hold on
plot(t,x(i:n*L+1:end)+y_ss);
end
hold off

figure % Pump voltage
hold on
plot(t,x(n*L+1:n*L+1:end)+u_ss);
hold off

%% Tradeoff curve between gamma and y

y = 3:0.5:6;

figure
for i = 1:length(y)
   d = [y(i); y(i); 5; 5; 8.18; 1.82];
   L = 10;
   [F, phi_Q, g, C_p, d_p, x, gam] = OptRobustControl(A,B,C,D,d,Q,R,S,T,L,N,sigma,u_0,x_0);
   u_rob = x(n*L+1:n*L+1:end);
   L = 1;
   [F, phi_Q, g, C_p, d_p, x, gam] = OptRobustControl(A,B,C,D,d,Q,R,S,T,L,N,sigma,u_0,x_0);
   u_nom = x(n*L+1:n*L+1:end);
   x_rob = x_0(1:2);
   x_nom = x_0(1:2);
   GammaRob = [];
   GammaNom = [];
   for j = 1:length(u_rob)
       x_rob = [x_rob A*x_rob(:,j)+B*u_rob(j)];
       x_nom = [x_nom A*x_nom(:,j)+B*u_nom(j)];
       GammaRob = [GammaRob (x_rob(:,j)+y_ss)'*Q(1:2,1:2)*(x_rob(:,j)+y_ss) + u_rob(j)'*R*u_rob(j) + 2*x_rob(:,j)'*S(1:2)*u_rob(j)];
       GammaNom = [GammaNom (x_nom(:,j)+y_ss)'*Q(1:2,1:2)*(x_nom(:,j)+y_ss) + u_nom(j)'*R*u_nom(j) + 2*x_nom(:,j)'*S(1:2)*u_nom(j)];
   end
   GammaRob = sum(GammaRob);
   GammaNom = sum(GammaNom);
   hold on
   plot(y(i),GammaRob,'b--o',y(i),GammaNom,'r--o');
end
