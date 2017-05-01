function [ vI,vJ, vI_dot, vJ_dot, T_cutoff ] = NetTideEquation( beta, theta, beta_prime, alpha, gamma, N, vT, I0, J0)
% Author: chengxi zang
% Date: Jan. 11, 2016
% Input: modeling parameters, time epocs vector, and initial cumulative I and J values.
% Output: the predicted I and J values.
%
% Goal: see whether we can have a model
% that will automatically create power law growth of node, link and densification.
% 
% Approach: 
% I-dot = beta * I / t^theta * (N - I )
% J-dot = beta_prime * I / t^theta ( alpha * (I-1)^gamma - J/I ) + 2*I-dot
% 
% where I(t) = # of infected nodes
%       J(t) = # of infected edges
%
% Assumption:
% '''
% 


lensOfT = length(vT);
vI = zeros(lensOfT,1);
vJ = zeros(lensOfT,1);
vI_dot = zeros(lensOfT,1);
vJ_dot = zeros(lensOfT,1);

tic = 1;
vI(tic) = I0;
vJ(tic) = J0;
tic = tic + 1;

I = I0;
J = J0;

while tic <= lensOfT
   
    t = vT(tic);
    I_dot =  beta * I / t^theta * (N - I);
    if alpha*(I-1)^gamma - J/I > 0
        J_dot =  2*beta_prime * I / t^theta *(alpha*(I-1)^gamma - J/I) + I_dot*2;
    else
        J_dot = I_dot*2;
    end

    I = I + I_dot;
    J = J + J_dot; 
    if (I > N)
        break;
    end
    if  I_dot < 10^(-4) 
        break
    end

    vI(tic) = I;
    vJ(tic) = J;
    vI_dot(tic) = I_dot;
    vJ_dot(tic) = J_dot;
    tic = tic+1;
end

T_cutoff = tic-1;

end

