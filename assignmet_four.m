% SC42135 SPECTRAL ANALYSIS OF NONLINEAR/INFINITE-DIMENSIONAL SYSTEMS
% Practice session 4 Dynamic mode decomposition
close all; clear all; clc

%% Task 1. Compute data matrices Y and Y'
m = 99;
ts = 2/m;
x0_1 = [1;0];
x0_2 = [0;1];

% Compute the states
[~,x1] = ode45(@P4T1_fx, (0:ts:ts*m), x0_1);
[~,x2] = ode45(@P4T1_fx, (0:ts:ts*m), x0_2);

% Domains for Y and Yp
dom1 = 1:length(x1)-1;
dom2 = 2:length(x1);

% Data amtrices
Y1 = [x1(dom1,1), x1(dom1,2)-x1(dom1,1).^2];
Y2 = [x2(dom1,1), x2(dom1,2)-x2(dom1,1).^2];
Y1p = [x1(dom2,1), x1(dom2,2)-x1(dom2,1).^2];
Y2p = [x2(dom2,1), x2(dom2,2)-x2(dom2,1).^2];

%% Task 2. Combine the matrices and compute A
Y = [Y1;Y2];
Yp = [Y1p;Y2p];

A = Yp'*pinv(Y');

%% Task 3. Compute Ac. NOT VERIFIED YET
Ac = logm(A)/ts;


%% Task 5. Find matrix C s.t. x[n]=Cy[n]
d=2;

% YY' data matrices using new measurement functions
Yn1 = P4T4(x1,d);
Yn2 = P4T4(x2,d);
Yn1p = Yn1(dom2,:);
Yn1 =  Yn1(dom1,:);
Yn2p = Yn2(dom2,:);
Yn2 =  Yn2(dom1,:);

Yn = [Yn1;Yn2];
Ynp = [Yn1p;Yn2p];

% New A matrix
An = Ynp'*pinv(Yn');

Xp = x1(dom2,:);
Xp2 = x2(dom2,:);
C = Xp'*pinv(An*Yn1');
%% Task 6. State trajectory for x[0] = [-0.7 1.1]

% ODE method
x0_t6 = [-0.7;1.1];
[~,xt6] = ode45(@P4T1_fx, (0:ts:ts*m), x0_t6);

% Koopman
y0 = P4T4(x0_t6,d);
x_hat = C*An*y0;

%% Task 7. Repeat tasks 5 & 6 with d=3
