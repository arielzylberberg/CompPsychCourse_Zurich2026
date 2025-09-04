addpath(genpath('../functions_addtopath/'))

load('../data/test_data.mat','rt','coh','c','choice');

D = struct('rt',rt,'coh', coh,'choice',choice,'c',c);

%% one example run
kappa = 15;   % signal-to-noise
ndt_m = 0.2;  % non-decision time, mean [s]
ndt_s = 0.01; % non-decision time, standard dev. [s]
B0    = 2;    % bound height
a     = 0.5;  % bound collapse parameter
d     = 0;    % bound parameter, time of decay [s]
y0    = 0;    % bias, on starting point
coh0  = 0;    % bias, on coherence

theta = [kappa,ndt_m,ndt_s,B0,a,d,coh0,y0];

pars = struct('plot_flag',true);

[err,P] = wrapper_DTB_parametricbound(theta,D,pars);

%% fitting

% kappa, ndt_mu, ndt_sigma, B0, a, d, coh0, y0
tl = [5,  0.1, .01 ,0.5  , -1, -3,0,0];
th = [40, 0.7, .15 ,4    , 4 ,4,0,0];
tg = [15, 0.2, .02 ,1    , 0.1 ,1,0,0];

pars = struct('plot_flag',true,'USfunc','Logistic');

fn_fit = @(theta) (wrapper_DTB_parametricbound(theta,D,pars));

MaxFunEvals = 100; % For the tutorial only, so it does not take too long
options = optimset('Display','final','TolFun',.01,'FunValCheck','on',...
    'MaxFunEvals',MaxFunEvals);


ptl = tl;
pth = th;
[theta, fval, exitflag, output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

% run the DTB with the best fit
[~,P] = fn_fit(theta);

% save fitted params
% save('fit_output','theta','fval','exitflag','output','tl','th','tg','pars');

%% evaluate best solution
load ../data/test_data.mat
load fit_output
D = struct('rt',rt,'coh', coh,'choice',choice,'c',c);
pars = struct('plot_flag',true,'USfunc','Logistic');
fn_fit = @(theta) (wrapper_DTB_parametricbound(theta,D,pars));
fn_fit(theta);



