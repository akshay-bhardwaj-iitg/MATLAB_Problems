close all
clear all

numSims = 5;    %5 5 runs
tBounds = [0 10];   % Time bound
N = 1000;
dt = (tBounds(2) - tBounds(1))/N;
b_init = 100;   % Initial Bt value

pd = makedist('Normal','mu',0,'sigma',sqrt(dt));

c = [0.7, 1.5, 0.06];   % K1, s1, sigma

ts = linspace(tBounds(1), tBounds(2), N); % From t0-->t1 with N points
bs = zeros(1,N); % 1xN Matrix of zeros

bs(1) = b_init;
%% Computing the Process
for j = 1:numSims
    for i = 2:numel(ts)
        t = tBounds(1) + (i-1).*dt;
        x = bs(i-1);
        a = -c(1).*x + c(2) + 0.5*c(3)*c(3).*x;
        b = -c(3).*x;
        dW = random(pd);

        bs(i) = x + a.*dt + b.*dW;
    end
       
    plot(ts, bs, 'r')
    hold on;
    xlabel('Time (s)');
    ylabel('B(t)');
    title('Solution estimated by Eulers method')
end
   
       
