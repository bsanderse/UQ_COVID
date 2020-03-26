function [QoI,solution] = solve_corona(X,P)
%
%   This code is can be run using
%        Matlab (commercial package)
%        Octave (open source version)
%
%   parameters for SEIR model
%
%   S = susceptible population
%   E = Exposed (infected, not yet infectious)
%   I = Infectious (now can infect others)
%   R = Removed (got sick, now recovered and immune, or died :( )
%   N = total population = (S + E + I + R)
%
%
%      note: added cRI/N term:  disease
%             mutates, can cause reinfection, or immunity lost
%             This assumes that mutated form jumps to Infected population
%             Can also assume that mutated form jumps to Exposed population
%             For now, we assume c=0 (no mutation has been observed)
%
%    dS/dt = -beta*S*I/N
%    dE/dt = +beta*S*I/N - sigma*E
%    dI/dt = +sigma*E -gamma*I + c*R*I/N
%    dR/dt = gamma*I -c*R*I/N

%   this file passes "seir.m" function to ode solver
%     ode systen is specified in "seir.m" file


%
%   Parameters from:
%
%   Wang, H., Wang, Z., Dong, Y. et al. Phase-adjusted estimation of
%   the number of Coronavirus Disease 2019 cases in Wuhan, China.
%   Cell Discov 6, 10 (2020). https://doi.org/10.1038/s41421-020-0148-0

%
%   Wuhan, Feb 2020
%          Based on estimates for original outbreak in Wuhan
%
%      These parameters are pretty much guestimates, but are probably
%      the right order of magnitude
%


% uncertainty vector X:

% initial conditions:
S_0 = X(1);
E_0 = X(2);
I_0 = X(3);
R_0 = X(4);
N   = S_0 + I_0 + E_0 + R_0;  % N = total population

% model parameters:
params.c     = X(5);
params.sigma = 1./X(6);
params.gamma = 1./X(7);
params.N     = N;

R_mult = X(8);

% other parameters in P:
t_start = P.t_start;
t_end   = P.t_end;

% R_zero = number of people infected by each infectious person
%          this has nothing to do with "R" = removed above
%          or R_0 (initial value of recovered)
%          but is common terminology (confusing, but usual notation)
%     time dependent, starts offf large, than drops with
%         time due to public health actions (i.e. quarantine, social distancing)
%
%    R_zero > 1, cases increase
%    R_zero < 1; cases peak and then drop off

%   R_zero declining with time https://www.nature.com/articles/s41421-020-0148-0

%   beta = R_zero*gammma (done in "seir.m" )


R_zero_array = zeros(6,2);
%  table of:   time(days)  R_zero
%            ....     ....
%            ....     ....
%            ....     ....
% linearly interpolate between times
%
% Note: this is different from Wang et al (2020), which assumes
%       piecewise constant values for R_zero
%

R_zero_array =        [0.0  3.0  ; ...     % t=0 days;    R_zero = 3.0
    60.0  2.6 ; ...      % t = 60 days; R_zero = 2.6
    70.0  1.9 ;...       % t = 70 days; R_zero = 1.9
    84.0  1.0; ...       % t = 84 days; R_zero = 1.0
    90.0  .50;...        % t = 90 days; R_zero = .50
    1000. .50 ] ;        % t = 1000 days; R_zero =.50

R_zero_array(:,2) = R_zero_array(:,2)*R_mult;
% R_zero_array = R_zero_array*R_mult;

params.R_zero_array = R_zero_array;

%
%
%  time units = days

%
%

% [start_time  end_time] (days)
tspan = [t_start t_end];  % time in days


%

%    y(1) = S
%    y(2) = E
%    y(3) = I
%    y(4) = R


yinit = zeros(4,1);
yinit(1) = S_0;
yinit(2) = E_0;
yinit(3) = I_0;
yinit(4) = R_0;



tol = 1.e-6;  % ode solver tolerance

%
%

% set 'Stats','on' to get more info

options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats','on');

%
%    note: set Refine switch to avoid interpolation
%
%options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats','on','Refine',1);

runtime = cputime ;

%    non-stiff solver, matlab and Octave
[t,y] = ode45(@(t, y) seir(t,y, params) , tspan, yinit, options);

%
%     stiff solver, matlab only
% [t,y] = ode15s( @(t, y) seir(t,y, params), tspan, yinit, options);
%

runtime =  cputime - runtime;

nsteps = length(t);

%disp(sprintf('number of steps:\t%15.5f',nsteps));
total_cases(:,1) = y(:,2) + y(:,3) + y(:,4);
total_cases_active(:,1) = y(:,2) + y(:,3) ;

if (P.plot_solution == 1)
    figure(11)
    hold on
    
    subplot(2,1,1), plot( t, y(:,1),'b-');
    xlabel('time(days)');
    ylabel('S: susceptible');
    
    hold on
    
    subplot(2,1,2), plot( t, y(:,2),'b-');
    xlabel('time(days)');
    ylabel('E: exposed');
    %
    figure(12)
    hold on
    subplot(2,1,1), plot( t, y(:,3),'b-');
    xlabel('time(days)');
    ylabel('I: infectious');
    
    hold on
    
    subplot(2,1,2), plot( t, y(:,4),'b-');
    xlabel('time(days)');
    ylabel('R: recovered');
    %
    
    figure(13)
    hold on
    
    plot( t, total_cases(:,1),'b-');
    xlabel('time(days)');
    ylabel('Total Cases: E+I+R ');
    %
    figure(14)
    hold on
    
    plot( t, total_cases_active(:,1),'b-');
    xlabel('time(days)');
    ylabel('Total Active Cases: E+I ');
    
end


S_end = y(nsteps, 1);
E_end = y(nsteps, 2);
I_end = y(nsteps, 3);
R_end = y(nsteps, 4);

total = S_end + E_end + I_end + R_end;

% QoI = y(end,:)';
QoI = max(total_cases_active);

solution.y = y;
solution.t = t;

% disp(sprintf('CPU time(sec):\t%15.5f',runtime));
%
% disp(sprintf('\n time (days): \t%10.2f \n', t(nsteps) ) );
%
% disp(sprintf('total population:\t%10.2f',total));
%
% disp(sprintf('initial infected:\t%10.2f',I_0));
%
% disp(sprintf('\n total cases (E+I+R) at t= %10.2f: %10.2f \n',...
%                  t(nsteps), E_end + I_end + R_end ));
%
%
% disp(sprintf('\n Recovered at t= %-10.2f: %10.2f \n', t(nsteps) ,R_end));
% disp(sprintf('\n Infected (infectious) at t= %-10.2f: %10.2f \n', t(nsteps),I_end));
% disp(sprintf('\n Exposed (non-infectious) at t= %-10.2f: %10.2f \n',t(nsteps),E_end));
% disp(sprintf('\n Susceptable at t= %-10.2f: %10.2f \n', t(nsteps), S_end));




