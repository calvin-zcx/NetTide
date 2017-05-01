function [ paras, simuResult ] = fitNetTide( tN, vN, vN_dot, tE, vE, vE_dot, threshold, ininitialPara)
%Learning the parameters of NetTide by minimizing the least square
%errors of n(t) and e(t).
%   
    %% Set the learning algorithm.
    options.Algorithm = 'levenberg-marquardt';
    
    %% Set the learning process starts from first threshold nodes.
%     threshold = 0.5*10^-2;
    tick0 = find(vN>=vN(end)*threshold, 1, 'first'); %start
    if tick0 ~= 1
        tick0 = tick0-1;
    end
    
    %% Learn NetTide-Node.
    if isempty(ininitialPara)
        beta0 = 0.001; 0.001; %beta
        theta0 = 1;     %theta
        N0 = 13000;    % population potential15625;
        
        beta_prime0 = 1;.5; 1;0.1;0.15; 0.15;
        alpha0 = 1;0.1;0.001;0.01;
        gamma0 =  1;.5;1;
    else
        beta0 = ininitialPara.beta;
        theta0 = ininitialPara.theta;
        N0 = ininitialPara.N;
        
        beta_prime0 = ininitialPara.beta_prime;
        alpha0 = ininitialPara.alpha;
        gamma0 =  ininitialPara.gamma;
    end
    
    [XI,resnormI] = lsqnonlin(@(x) Objective_func_NetTideNode( x, tN(tick0:end), vN(tick0:end)), ...
        [beta0, theta0, N0],[0,0,0],[1, Inf, Inf],options);        
    beta_learned = XI(1);
    theta_learned = XI(2);
    N_learned = XI(3);
 
 
%%    Learn NetTide-Link.  
    [XJ, resnormJ] = lsqnonlin(@(x) Objective_func_NetTideLink( x, theta_learned,  tN(tick0:end), vN_dot(tick0:end), vN(tick0:end), vE(tick0:end)), ...
    [ beta_prime0, alpha0, gamma0],[0, 0, 0],[1, 1, 1],options);
         
    beta_prime_learned = XJ(1);
    alpha_learned = XJ(2);
    gamma_learned = XJ(3);

    %% Simulation results.
    
    I0 = vN(tick0);
    J0 = vE(tick0);
    [ vN_simu, vE_simu, vN_dot_simu, vE_dot_simu, T_cutoff ] =...
            NetTideEquation(  beta_learned, theta_learned, ...
                        beta_prime_learned, alpha_learned, gamma_learned, ...
                        N_learned, tN(tick0:end), I0, J0);
    %% Return results.                
    paras.beta_learned = beta_learned;
    paras.theta_learned = theta_learned;
    paras.N_learned = N_learned;
    paras.beta_prime_learned = beta_prime_learned;
    paras.alpha_learned = alpha_learned;
    paras.gamma_learned = gamma_learned;
    
    simuResult.tick0 = tick0;
    simuResult.vN_simu = vN_simu;
    simuResult.vE_simu = vE_simu;
    simuResult.vN_dot_simu = vN_dot_simu;
    simuResult.vE_dot_simu = vE_dot_simu;
    simuResult.T_cutoff = T_cutoff;
    simuResult.resnormI = resnormI;
    simuResult.resnormJ = resnormJ;
    
end
