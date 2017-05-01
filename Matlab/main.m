% Author: chengxi zang
% Date: May. 1, 2017
%
% Goal: Load the growth dynamics generated from two generators, 
% fit, plot and simulate  them.
%
% 

close all
clear

%% Read in the data
%% Choose the data generated from two generators
nodeFileName = 'nettide_survival\\nodeAdd_rate';
linkFileName = 'nettide_survival\\linkAdd_rate';

% nodeFileName = 'nettide_process\\nodeAdd_rate';
% linkFileName = 'nettide_process\\linkAdd_rate';

Data = readin_data(nodeFileName, linkFileName);

%% Set the color order for each selected instance for plottNng
% mcolor = get(gca,'colororder');
mcolor = [0,0.447,0.741;0.850,0.325,0.098;0.929,0.694,0.125;0.494,0.184,0.556;0.466,0.674,0.188;0.301,0.745,0.933;0.635,0.078,0.184];
vparas = [];
vexponent = [];
 
f = figure;
%for dynamic simulatNon [1:6 9:11 15]
%for dynamic simulatNon [2 3 5:8 14 17 18 20 ]
iter = 0;

%% Selected instances for illustration
%%NettNde-SurvNval:[1:3 5 10:13 17 22]
%%NettNde-Process:[10 12 14 1 3 5:8 15 ]
for dim =  [1:3 5 10:13 17 22];
% for dim = [10 12 14 1 3 5:8 15 ];
    iter = iter + 1; 
   %% Prepare data 
    tN = Data{dim,1}; %time of n(t)
    vN_dot = Data{dim,2};  %node rate
    tE = Data{dim,3}; % time of e(t)
    vE_dot = Data{dim,4}; %edge rate
    vN = cumsum(vN_dot); %n(t)
    vE = cumsum(vE_dot); %e(t)
    % Log transformed data
    log_vN_dot = log(vN_dot);
    log_vE_dot = log(vE_dot);
    log_vN = log(vN);
    log_vE = log(vE);

 %% Learn for initNal parameters
    [ paras, sr ] = fitNetTide( tN, vN, vN_dot, tE, vE, vE_dot,  10^-2, []);

    ininitialPara.beta = paras.beta_learned;
    ininitialPara.theta= paras.theta_learned;
    ininitialPara.N = paras.N_learned;
    ininitialPara.beta_prime = paras.beta_prime_learned;
    ininitialPara.alpha = paras.alpha_learned;
    ininitialPara.gamma = paras.gamma_learned;

%% Learn for final result
    [ paras, sr ] = fitNetTide( tN, vN, vN_dot, tE, vE, vE_dot,  10^-2, ininitialPara);
    tick0 = sr.tick0;
    vN_simu = sr.vN_simu;
    vE_simu = sr.vE_simu;
    vN_dot_simu = sr.vN_dot_simu;
    vE_dot_simu = sr.vE_dot_simu;
    T_cutoff = sr.T_cutoff;
    resnormI = sr.resnormI;
    resnormJ = sr.resnormJ;

    beta_learned = paras.beta_learned;
    theta_learned = paras.theta_learned;
    N_learned = paras.N_learned;
    beta_prime_learned = paras.beta_prime_learned;
    alpha_learned = paras.alpha_learned;
    gamma_learned = paras.gamma_learned;


%% keep all the learned parameters for further analysis
    %1 beta_learned, 2theta_learned, 3N_learned, 
    %4beta_prime_learned, 5alpha_learned, 6gamma_learned, 7vN(tick0), 8vE(tick0)];
    vparas = [vparas; beta_learned, theta_learned, N_learned, beta_prime_learned, alpha_learned, gamma_learned, vN(tick0), vE(tick0), resnormI, resnormJ];
    ix = mod(iter, size(mcolor,1))+1;

%% Fit the power law early growth
    a = tN(10:100);%(tick0:end);
    b = vN(10:100);%(tick0:end);
    loga =  log(a);
    logb = log(b);
    paraI = polyfit(loga, logb, 1);

    a = tE(10:100);
    b = vE(10:100);
    loga =  log(a);
    logb = log(b);
    paraJ = polyfit(loga, logb, 1);
    vexponent = [vexponent; paraI(1), paraJ(1) ];
    
%% Node cumulative Plot Demo
    plot(tN(1:10:end), vN(1:10:end), 's', 'color',  mcolor(ix,:), 'markersize', 5);%'color',[.35 .35 .35] 
    hold on
    plot(tN(tick0:end), vN_simu, 'color', mcolor(ix,:),   'linewidth', 1.5);%, 'color', mcolor(ix,:)
    hold on
%
%% Node rate Plot Demo
%     plot(tN(1:5:end), vN_dot(1:5:end), '-.', 'color', mcolor(ix,:));%, 'color', mcolor(ix,:), 'color',[.5 .5 .5]);%, 'markersize', 10
%     hold on
%     plot( tN(tick0+1:end), vN_dot_simu(2:end), 'color', mcolor(ix,:), 'linewidth', 1.5);%, 'color', mcolor(ix,:)'r', 'color', mcolor(dim,:),
%     hold on;
%     
%% Link cumulative Plot Demo.
%     plot(tE(1:10:end), vE(1:10:end), 'o', 'color', mcolor(ix,:), 'markersize', 5);%[.35 .35 .35] 
%     hold on
%     plot( tE(tick0:end), vE_simu, 'color', mcolor(ix,:),  'linewidth', 1.5);%, 'color', mcolor(ix,:)
%     hold on;
%
%% Link rate Plot Demo
%     plot(tE(1:2:end), vE_dot(1:2:end), '-.', 'color', mcolor(ix,:));%, 'color', mcolor(ix,:), 'color',[.5 .5 .5]);%, 'markersize', 10
%     hold on
%     plot( tE(tick0+1:end), vE_dot_simu(2:end), 'color', mcolor(ix,:), 'linewidth', 1.5);%, 'color', mcolor(ix,:)'r', 'color', mcolor(dim,:),
%     hold on;
%
%% DensificatNon
%     loglog(vN(1:end), vE(1:end), 'o', 'color',[.35 .35 .35] , 'markersize', 5);
%     hold on;
%     loglog( vN_simu, vE_simu, 'color', mcolor(ix,:),  'linewidth', 1.5);%, 'color', mcolor(ix,:)
%     hold on;
%    
end

grid
set(gca,'FontSize',13);
xlabel('time', 'fontsize', 24)

%% Choose the corresponding label and legend
ylabel('n(t)', 'fontsize', 24)
% ylabel('dn(t)/dt', 'fontsize', 24)
% ylabel('e(t)', 'fontsize', 24)
% ylabel('de(t)/dt', 'fontsize', 24)

legend({'NT-Surv', 'NetTide'}, 'fontsize', 20);
% legend({'NT-Proc', 'NetTide'}, 'fontsize', 20);


set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% fname = 'nettide_survival//figure//simuNt' ;
% fname = 'nettide_survival//figure//simuNtdot' ;
% fname = 'nettide_survival//figure//simuEt' ;
% fname = 'nettide_survival//figure//simuEtdot' ;
% print(f, fname, '-dpdf','-r0')
% print(f, fname, '-djpeg','-r0')