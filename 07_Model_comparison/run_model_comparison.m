addpath(genpath('../functions_addtopath/'));

%% load data
D = load('../data/roitman_data.mat');
fit_nobias = load('../05_Fitting_exercise_roitman/optim.mat');
fit_bias = load('../05b_Fit_Roitman_data_with_bias/optim.mat');

%% plot the model without bias
% simulate on a finer grid
coh_fine = linspace(-0.6,0.6,101); % finer coherences, for nicer plots
Dfine = struct('coh',coh_fine);
params_fine = struct('plot_flag',1,'USfunc','Logistic');
[~, Pfine] = wrapper_DTB_parametricbound(fit_nobias.theta, Dfine, params_fine);

% calc non decision time
ndt_m = fit_nobias.theta(2);
ndt_s = fit_nobias.theta(3);
pd = makedist('Normal','mu',ndt_m,'sigma',ndt_s);
pd_trunc = truncate(pd,0,inf);
ndt_mean = pd_trunc.pdf(Pfine.t)*Pfine.t(:)/sum(pd_trunc.pdf(Pfine.t));

p = plot_finer_coh_grid(D, coh_fine, Pfine, ndt_mean);
% p.current_ax(1);
% ht = title('Model with no bias term');
% set(ht,'FontWeight','normal');

p.unlabel_center_plots();
p.format('FontSize',20);

%% plot the model with bias
% simulate on a finer grid

[~, Pfine] = wrapper_DTB_parametricbound(fit_bias.theta, Dfine, params_fine);

% calc non decision time
ndt_m = fit_bias.theta(2);
ndt_s = fit_bias.theta(3);
pd = makedist('Normal','mu',ndt_m,'sigma',ndt_s);
pd_trunc = truncate(pd,0,inf);
ndt_mean = pd_trunc.pdf(Pfine.t)*Pfine.t(:)/sum(pd_trunc.pdf(Pfine.t));

p = plot_finer_coh_grid(D, coh_fine, Pfine, ndt_mean);
% p.current_ax(1);
% ht = title('Model with bias term');
% set(ht,'FontWeight','normal');

p.unlabel_center_plots();
p.format('FontSize',20);

%% AIC and BIC
% aic = -2*logl + 2*numParams;

numParams_nobias = 6;
numParams_bias = 7;

aic_nobias = 2*fit_nobias.fval + 2*numParams_nobias;
aic_bias = 2*fit_bias.fval + 2*numParams_bias;

% bic = -2*logl + numParams*log(numObs);
numObs = sum(~isnan(D.rt));
bic_nobias = 2*fit_nobias.fval + numParams_nobias*log(numObs);
bic_bias = 2*fit_bias.fval + numParams_bias*log(numObs);



