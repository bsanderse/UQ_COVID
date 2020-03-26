% Wuhan parameters

%% initial conditions
S_0 = 11.0e+6;  % Wuhan city  excluding initial infected, exposed population, 

I_0 = 40.0;   % initial infected from market
                  %https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news–wuhan-coronavirus/
                  
E_0 = 20.*I_0;    % initial exposed
                      % https://www.medrxiv.org/content/10.1101/2020.01.23.20018549v1.full.pdf
                  
R_0 = 0;         % initial recovered (not to be confused with R_zero, below)
                     % initially, no one has recovered

%% time integration
t_start = 0;
t_end   = 180; % in days
                     
%% model parameters

c = 0.0;     % no mutation (yet)
             % maybe this happens later?

sigma = 1./5.2;   %  https://doi.org/10.1056/NEJMoa2001316 (2020).

gamma = 1./18.; % https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/news–wuhan-coronavirus



%% R_zero array

% R_zero_array =    [0.0  3.0  ; ...     % t=0 days;    R_zero = 3.0
%                    60.0  2.6 ; ...      % t = 60 days; R_zero = 2.6
%                    70.0  1.9 ;...       % t = 70 days; R_zero = 1.9
%                    84.0  1.0; ...       % t = 84 days; R_zero = 1.0
%                    90.0  .50;...        % t = 90 days; R_zero = .50
%                    1000. .50 ] ;        % t = 1000 days; R_zero =.50


%% define uncertainties in X
% X(1) = S_0;
% X(2) = E_0;
% X(3) = I_0;
% X(4) = R_0;
% X(5) = c;
% X(6) = sigma or 1/sigma
% X(7) = gamma or 1/gamma
% X(8) = Rzero multiplier

ndim = 8;

% marginal distribution X1
Input.Marginals(1).Name = 'S_0';
Input.Marginals(1).Type = 'Constant'; 
Input.Marginals(1).Parameters = S_0;

% marginal distribution X2
Input.Marginals(2).Name = 'E_0';
Input.Marginals(2).Type = 'Constant'; 
Input.Marginals(2).Parameters = E_0; %[0.9*E_0 1.1*E_0]; 
% Input.Marginals(2).Bounds = [0.9*E_0 1.1*E_0]; 

% marginal distribution X3
Input.Marginals(3).Name = 'I_0';
Input.Marginals(3).Type = 'Constant'; 
Input.Marginals(3).Parameters = I_0; %[0.9*I_0 1.1*I_0]; % scale and shape parameter
% Input.Marginals(3).Bounds = [0.9*I_0 1.1*I_0]; % scale and shape parameter

% marginal distribution X4
Input.Marginals(4).Name = 'R_0';
Input.Marginals(4).Type = 'Constant'; 
Input.Marginals(4).Parameters = R_0;

% marginal distribution X5
Input.Marginals(5).Name = 'c';
Input.Marginals(5).Type = 'Constant'; 
Input.Marginals(5).Parameters = c;

% marginal distribution X6
Input.Marginals(6).Name = 'sigma_{inv}';
Input.Marginals(6).Type = 'Uniform'; 
Input.Marginals(6).Parameters = [4 7];

% marginal distribution X7
Input.Marginals(7).Name = 'gamma_{inv}';
Input.Marginals(7).Type = 'Uniform'; 
Input.Marginals(7).Parameters = [15 22];

% marginal distribution X8
Input.Marginals(8).Name = 'Rzero';
Input.Marginals(8).Type = 'Uniform'; 
Input.Marginals(8).Parameters = [0.9 1.1];
