function  yprime = seir(t,y, params)

%
%    dS/dt = -beta*S*I/N
%    dE/dt = +beta*S*I/N - sigma*E 
%    dI/dt = +sigma*E -gamma*I + c*R*I/N
%    dR/dt = gamma*I -c*R*I/N
%
%
%
%   yprime = [dS/dt  dE/dt dI/dt   dRdt]'
%
%   input:
%       t  current time
%       y vector of current soln values
%         y(1) = S, y(2) = E, y(3) = I, y(4) = R
%
%     parameters in "params" 
%           beta, N, sigma, gamma, c, R_zero_array (table of values) 
%
%   output: (col vector)
%    yprime(1) = dS/dt
%    yprime(2) = dE/dt
%    yprime(3) = dI/dt
%    yprime(4) = dR/dt


   R_zero_array = params.R_zero_array;
           % Table of: time    R_zero
           %            ...     ....
           %            ...     ....


    min_t = R_zero_array(1,1);
      n_table = length( R_zero_array(:,1) );
    max_t = R_zero_array(n_table,1);
     t_val = max( min_t, min( t, max_t) );

   R_zero = interp1( R_zero_array(:,1), R_zero_array(:,2), t_val);

   gamma = params.gamma;

   beta = R_zero*gamma; 
 
   N = params.N;
   sigma = params.sigma;
   c = params.c;
%
%
     S = y(1);
     E = y(2);
     I = y(3);
     R = y(4);

     yprime = zeros(4,1);

     yprime(1) = -beta*S*I/N;
     yprime(2) = +beta*S*I/N - sigma*E;
     yprime(3) = +sigma*E -gamma*I + c*R*I/N;
     yprime(4) = gamma*I -c*R*I/N;
