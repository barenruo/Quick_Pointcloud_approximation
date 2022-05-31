clc
clear all;
close all;
%----------------Simulate data-------------------------%
% functional model
%f=@(x) cos(20*x) + exp(1.5*x);
f=@(x) cos(50*x) + exp(1.5*x);
%initial value of polynomial degree
p=80;
%vector for plotting the original function
x_print = linspace(-1,1,1000);
x_print = x_print';
y_print = f(x_print);
%vector for evaluating the approxiamation

 figure(2)
 plot(x_print,y_print,'Markersize',5)

hold on

%adjusting the polynomial degree
%initialise parameters
accuracy=10^-17;
coef=inf;
%while min(abs(coef))> accuracy
sigma = zeros(14,1);
last = inf;
while last > accuracy
j = 0:p;
xfft = cos(j*pi/p)';%the chebyshev points
yfft = f(xfft);     %the associate y values
udv = [yfft;yfft(end-1:-1:2)]; %fill the rest area on the unit disk
fc = real(fft(udv)); %fft
coef = (fc(1:p+1))/p;
coef(1) = coef(1)/2;
coef(end) = coef(end)/2;%calculate the coefficients
last = abs(coef(end));
p=p+1;
end
yp_d = zeros(size(x_print));
for i = 0:p-1
   yp_d = yp_d + coef(i+1)*cos(i*acos(x_print));
end
  num_pts=p; %number of points 
x=linspace(0,pi,num_pts)';
x = -cos(x);%chebyshev points
y = f(x); 
figure(2)
plot(x,y,'.','Markersize',5)
hold on
figure(2)
plot(x_print,yp_d,'Markersize',5)
legend('Function value','Control points','Aproximation')
