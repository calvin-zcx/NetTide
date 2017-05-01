function [ obj ] = Objective_func_NetTideNode( X, vT, vI_empirical)
% Obejective function for the NetTide-Node
% Author: Chengxi Zang
% Date: Jan. 11, 2016
% 
% % I-dot = beta /t^theta* I * (N - I )
% % A discrete time version


beta = X(1);
theta = X(2);
N = X(3);

lensOfT = length(vT);
tic = 1;
I0 = vI_empirical(tic);

vI = zeros(lensOfT,1);
vI_dot = zeros(lensOfT,1);

vI(tic) = I0;
tic = tic+1;
I = I0;
tend = lensOfT;
while tic <= lensOfT
%     gammaS = 1;
%     I_dot = beta * I^gammaI * ( N - I )^1;
%     I_dot = beta * I^1 * ( N - I )^1;
    t = vT(tic);
    I_dot = beta * I / t^theta * ( N - I );
    I = I + I_dot;

    if (I > N)
%         tend = tic;
        break;
    end

    vI_dot(tic) = I_dot;
    vI(tic) = I;

    tic = tic + 1;
end
    x = vI(1:tend);
    y = vI_empirical(1:tend);
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




