% call corona solver
clear all
close all
clc

input_file = 'Wuhan.m';

run(input_file);

%% uncertainties in X
X(1) = S_0;
X(2) = E_0;
X(3) = I_0;
X(4) = R_0;
X(5) = c;
X(6) = sigma;
X(7) = gamma;


%% parameters in P
P.t_start = t_start;
P.t_end = t_end;



%% run model
[QoI,solution] = solve_corona(X,P);

figure(1)
plot(solution.t,solution.y(:,2) + solution.y(:,3))