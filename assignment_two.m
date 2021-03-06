% SC42135 SPECTRAL ANALYSIS OF NONLINEAR/INFINITE-DIMENSIONAL SYSTEMS
% Practice session 2
clc; clear all; close all;
%% Task 1
% polynomial coefficients
a_coeffs = [-2 3 1 2];
L = 15;
% matrix T
deriv_operator = zeros(L*2+1,1);
for ii = 1:size(deriv_operator,1)
    deriv_operator(ii) = (L+1-ii)*1i;
end
spectrumT  = zeros(size(deriv_operator));
for kk=1:size(a_coeffs,2)
    spectrumT = spectrumT + a_coeffs(kk)*deriv_operator.^(kk-1);
end
%% Task 2
d = [-pi,pi];
A = chebop(d);
A.bc = 'periodic';
A.op = @(u) a_coeffs(1)*u+a_coeffs(2)*diff(u)+a_coeffs(3)*diff(u,2)...
          + a_coeffs(4)*diff(u,3);
[vv,D] = eigs(A,5);
%% Task 3
x = chebfun('x', d);
u0 = sin(x)^2+0.5*cos(x);
t = [0.1 0.5 2];
u = expm(A, t, u0);
%% Task 4 [Repeat the whole thing with an(1)=2]
b_coeffs = [2 3 1 2];
B = A;
B.op = @(v) b_coeffs(1)*v+b_coeffs(2)*diff(v)+b_coeffs(3)*diff(v,2)...
          + b_coeffs(4)*diff(v,3);
[v1,D1] = eigs(B,5);      
v = expm(B, t, u0);
%% Extra task


%% PLOTS

figure
set(gcf,'color','w');
plot(spectrumT,'o')
xlabel('Real Axis')
ylabel('Imaginary Axis')
grid('on')
axis([-7 0 -15 15])

figure;
set(gcf,'color','w');
plot(diag(D), '*','MarkerSize', 6,'LineWidth',2)
xlabel('Real Axis')
ylabel('Imaginary Axis')
title('Linear Operator eigenvalues')
% axis([-50 0 -20 20])
grid('on')

figure;
title('Solutions at time t')
set(gcf,'color','w');
subplot(3,1,1)
hold on
plot(u{1})
plot(v{1})
title(sprintf('time = %f',t(1)))
subplot(3,1,2)
plot(u{2})
hold on
plot(v{2})
title(sprintf('time = %f',t(2)))
subplot(3,1,3)
plot(u{3})
hold on
plot(v{3})
title(sprintf('time = %f',t(3)))

figure
title('Solutions at time 2.0')
yyaxis right
p1 = plot(u{3});
p1.LineWidth = 2
yyaxis left
plot(v{3})
l = legend(sprintf('coeffs [%d %d %d %d]',a_coeffs), ...
           sprintf('coeffs [%d %d %d %d]',b_coeffs) );
l.Position = [0.7372    0.2302    0.0698    0.0379];
l.FontSize = 12;


figure;
title('Solutions at time t')
set(gcf,'color','w');
for ii = 1:6
    subplot(3,2,ii)
    if(mod(ii,2))
        plot(u{ceil(ii/2)})
    else
    plot(v{ceil(ii/2)})
    end
    title(sprintf('time = %f',t(ceil(ii/2))))
end

subplot(3,1,2)
plot(u{2})
hold on
plot(v{2})
title(sprintf('time = %f',t(2)))
subplot(3,1,3)
plot(u{3})
hold on
plot(v{3})
title(sprintf('time = %f',t(3)))

figure;
set(gcf,'color','w');
subplot(2,1,1)
plot(real(D), imag(D), '*','MarkerSize', 6,'LineWidth',2)
xlabel('Real Axis')
ylabel('Imaginary Axis')
title(sprintf('Linear Operator eigenvalues [%d %d %d %d]',a_coeffs))
axis([-7 2 -15 15])
grid('on')
subplot(2,1,2)
plot(real(D1), imag(D1), '*','MarkerSize', 6,'LineWidth',2)
xlabel('Real Axis')
ylabel('Imaginary Axis')
title(sprintf('Linear Operator eigenvalues [%d %d %d %d]',b_coeffs))
axis([-7 2 -15 15])
grid('on')

figure;
plot(u{1})
hold on
plot(u{2})
plot(u{3})
