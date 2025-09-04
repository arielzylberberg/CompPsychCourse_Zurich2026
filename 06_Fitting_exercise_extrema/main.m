addpath(genpath('../functions_addtopath/'))

load ../data/sim_dat_extrema.mat

%% fitting

% kappa, ndt_mu, ndt_sigma, B0, a, d, coh0, y0
tl = [5,  0.1, .01 ,0.5  , -1, -3,0,0];
th = [40, 0.7, .15 ,4    , 4 ,4,0,0];
tg = [15, 0.2, .02 ,1    , 0.1 ,1,0,0];

pars = struct('plot_flag',true,'USfunc','Logistic');

for i = 1:length(sim)

    fn_fit = @(theta) (wrapper_DTB_parametricbound(theta,sim(i),pars));

    MaxFunEvals = 50; % For the tutorial only, so it does not take too long
    options = optimset('Display','final','TolFun',.01,'FunValCheck','on',...
        'MaxFunEvals',MaxFunEvals);


    ptl = tl;
    pth = th;
    [theta, fval, exitflag, output] = bads(@(theta) fn_fit(theta),tg,tl,th,ptl,pth,options);

    % run the DTB with the best fit
    [~,P] = fn_fit(theta);
    
    theta_best(i,:) = theta;
end

% save optim theta_best MaxFunEvals

%% plot best fits

load optim
load ../data/sim_dat_extrema.mat
pars = struct('plot_flag',true,'USfunc','Logistic');
[~,P1] = wrapper_DTB_parametricbound(theta_best(1,:),sim(1),pars);
[~,P2] = wrapper_DTB_parametricbound(theta_best(2,:),sim(2),pars);



