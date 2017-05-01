function [ obj ] = Objective_func_NetTideLink( X, theta, vT, vI_dot_empirical, vI_empirical, vJ_empirical)
% Obejective function for the NetTide-Link
% Author: Chengxi Zang
% Date: Jan. 10, 2016
% 
% % J_dot = beta_prime/t^theta * I *(alpha*(I-1)^gamma - J/I) + I_dot*2;
% % A discrete time version


beta_prime = X(1);
alpha = X(2);
gamma = X(3);

lensOfT = length(vT);
tic = 1;
vJ = zeros(lensOfT,1);
vJ_dot = zeros(lensOfT,1);
J = vJ_empirical(1);

vJ(1) = J;

tic = tic+1;

while tic <= lensOfT
    
    t = vT(tic);
    I_dot = vI_dot_empirical(tic); 
    I = vI_empirical(tic);
    
%     J_dot = beta_prime * I / t^theta *(alpha*(I-1)^gamma - J/I) + I_dot*2;
    if alpha*(I-1)^gamma - J/I > 0
        J_dot =  2*beta_prime * I / t^theta *(alpha*(I-1)^gamma - J/I) + I_dot*2;
    else
        J_dot = I_dot*2;
    end

    J = J + J_dot; 
%     if (I > N)
%         break;
%     end

    vJ_dot(tic) = J_dot;
    vJ(tic) = J;
   
    tic =tic + 1;
end
    x = vJ(1:lensOfT);
    y = vJ_empirical(1:lensOfT);
    obj = LogL2Distance(x, y);%x-y; %
end

% function [val] = L2Distance(x,y)
%     val = sum((x-y).*(x-y));
% end

function [val] = LogL2Distance(x,y)
    a = zeros(size(x));
    a(x>1) = log(x(x>1));
    b = zeros(size(y));
    b(y>1) = log(y(y>1));
%     val = sum((a-b).*(a-b));
    val = a-b;
end




