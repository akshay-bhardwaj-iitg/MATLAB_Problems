%% Initial data and calculations
Ww = 9.81;
W = @(Gs,e) Ww*(Gs+0.2*(e))/(1+e);
Wsat = @(Gs,e) Ww*(Gs+e)/(1+e);
h = @(H,U) (H*U);

P = @(W,Wsat,h,H,thita,phi) ((W*(H-h)+(h*(Wsat-Ww)))*(cos(thita*pi/180))*(tan (phi*pi/180)))/(((W*(H-h))+(h*Wsat))*sin(thita*pi/180));

mu_phi  = 35;   
var_phi = .08;
muphi    = log((mu_phi^2)/sqrt(var_phi+mu_phi^2));
sigmaphi = sqrt(log(var_phi/(mu_phi^2)+1));

mu_thita  = 20;   
var_thita = .05;
muthita    = log((mu_thita^2)/sqrt(var_thita+mu_thita^2));
sigmathita = sqrt(log(var_thita/(mu_thita^2)+1));

a    = 1;           % threshold level
NSIM = 1000;
b    = zeros(1,NSIM);
c    = zeros(1,NSIM);
d    = zeros(1,NSIM); 
x    = zeros(1,NSIM);


%% MCS using normal computing
fprintf('MONTE CARLO SIMULATION : \n');
tic;
for i = 1:NSIM
    H = unifrnd(2,8);
    U = unifrnd(0,1);
  phi    = lognrnd(muphi,sigmaphi);
  thita    = lognrnd(muthita,sigmathita);
  Gs = unifrnd(2.5,2.7);
  e = unifrnd(0.3,0.6);
  b(i) = W(Gs,e);
  c(i) = Wsat(Gs,e);
  x(i) = h(H,U);
  d(i) = P(W(Gs,e),Wsat(Gs,e),h(H,U),H,thita,phi);
end
toc;   t1 = toc;
[ff1,xx1] = ecdf(d);              % estimate empirical CDF
pf        = mean(d<=a);           % failure probability
var_MCS   = pf*(1-pf)/NSIM;       % variance 
std_MCS   = sqrt(var_MCS);
fprintf('Failure probability: %7.8f +- %g \n\n', pf, std_MCS);