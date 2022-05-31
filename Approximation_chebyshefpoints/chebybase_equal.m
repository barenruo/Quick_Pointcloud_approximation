clc
clear all;
close all;
% function
f=@(x) cos(25*x) + exp(1.5*x);
%initial value of polynomial degree
p=20;
%vector for plotting the original function
x_print = linspace(-1,1,1000);
y_print = f(x_print);
%vector for evaluating the approxiamation
x_d = linspace(-1,1,1000000);
x_d = x_d';
y_d = f(x_d);
figure(2)
plot(x_print,y_print,'g','linewidth',4)
title({'Function approximation with equally distributed'...
    ;'control points using Chebyshev basis'})
xlabel('x(m)')
ylabel('y(m)')
hold on

%adjusting the polynomial degree
%initialise parameters
accuracy=5;
v=inf;
int_v = inf;
tic
while int_v> accuracy
    num_pts=p+1; %number of points 
x = linspace(-1,1,num_pts)';%select p+1 point
y = f(x);                  %associated function value
yp=zeros(size(x));
A = [];                  %Vandermonde matrix
for i=0:p
 A=[A cos(i*acos(x))]; 
end
X_hat=linsolve(A,y);
yp_d = zeros(size(x_d));
for i = 0:p
   yp_d = yp_d + X_hat(i+1)*cos(i*acos(x_d));
  end
v1=y_d-yp_d;
  int_v = integralresidual(x_d,v1);
  p=p+1;
end
toc
 figure(2)
plot(x,y,'k.','Markersize',15)
hold on
figure(2)
plot(x_d,yp_d,'r','linewidth',1)
legend('Function value','Selected points','Aproximation')