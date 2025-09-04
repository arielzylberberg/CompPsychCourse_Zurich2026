%% add functions to path
addpath(genpath('../functions_addtopath/'));

%% figure defaults
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 1);

%% load data
load ../data/test_data.mat

%% run one set of parameters
kappa = 15;         % signal-to-noise parameters (\mu = kappa x coh)
ndt_m = 0.2;        % non-decision time, mean [s]
ndt_s = 0;          % non-decision time, standard deviation [s]
B0    = 1.2;        % bound height
coh0  = 0;          % bias (in units of coherence)
y0    = 0;          % starting-point bias, as a fraction of bound height
ndt_m_delta = 0;    % difference in non-decision time between right and left choices
theta = [kappa,ndt_m,ndt_s,B0,coh0,y0,ndt_m_delta]; % model parameters

params = struct('plot_flag',1);

% data struct
D = struct('rt',rt,...          % response time
           'coh', coh,...       % motion coherence
           'choice',choice,...  % choice [0,1]
           'c',c);              % correct [0,1]

% eval the model parameters and plot
[~,P] = wrapper_analytic_DDM(theta,D,params);

%% fit parameters, sum of squared errors

% kappa, ndt_mu, ndt_sigma, B0, coh0, y0, ndt_m_delta
tl = [5,  0.1, 0 ,.5 ,0,0,0]; % lower limit
th = [40, 0.5, 0 ,3  ,0,0,0]; % upper limit
tg = [15, 0.2, 0 ,1  ,0,0,0]; % guess

params = struct('plot_flag',1,'optim_method',2); % sum of squared errors of RT

fn_fit = @(theta) (wrapper_analytic_DDM(theta,D,params));

MaxFunEvals = 200; % For the tutorial only, so it does not take too long
options = optimset('Display','final','TolFun',.01,'FunValCheck','on',...
    'MaxFunEvals',MaxFunEvals);
ptl = tl;
pth = th;
[theta_SE,fval,exitflag,output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

% save fitted params
% save('fit_output_SE','theta_SE','fval','exitflag','output','tl','th','tg','params');




