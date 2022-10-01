%% the realizations of 𝐵𝑡 using Eular-Maruyama scheme
%{
---------------------------------------------------------------------------
*Created by:                Date:            
 Akshay Kumar Bhardwaj      November 2021    
*Mail: 
 akshay_kumar@iitg.ac.in
*Istitute:
 Indian Institue of Technology, Guwahati
---------------------------------------------------------------------------
%}
A=10; 
t=0; 
fs=100;
omega=2*pi*2;

%%
theta=rand(1,10000)*2*pi;
x1=A*sin(omega*t+theta);
Rxx1=[]; 
tp=10;

for tau=-tp:1/fs:tp
tmp=A*sin(omega*(t+tau)+theta);
tmp=mean(x1.*tmp);
Rxx1=[Rxx1 tmp];
end

tau=-tp:1/fs:tp;
Rxx=A^2/2*cos(omega*tau);

%% 
t=0:1/fs:20-1/fs;
x2=A*sin(omega*t);
[Rxx2, tau2]=xcorr(x2,x2,tp*fs,'unbiased');
tau2=tau2/fs;

%%

subplot(2,1,1);
plot(tau,Rxx1,tau,Rxx,'linewidth',1)
xlabel('tau')
ylabel('Autocorrelation')
title('Ensemble Average')
grid on

subplot(2,1,2);
plot(tau2,Rxx2,tau,Rxx,'g','linewidth',1)
xlabel('tau')
ylabel('Autocorrelation')
title('Time Average')
grid on