function [err,P] = wrapper_analytic_DDM(theta, D, params)
% function [err,P] = wrapper_analytic_DDM(theta, D, params)
% fits the drift-diffusion model with flat bounds to choice and/or reaction
% time data.
%
% Input parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%     theta: vector of parameters:
%       kappa,
%       mean of non-decision time,
%       standard deviation of non-decision times,
%       bound height,
%       bias (offset to the motion coherence),
%       bias (offset to the starting point; as a proportion of bound height),
%       ∆non-decision time (extra time added to the non-dec time of rightward choices)
%
%     D: struct containing the experimental data: 
%       D.choice: 0 and 1 for left and right responses respectively [ntrials]
%       D.rt: response time for each trial [ntrials]
%       D.coh: signed motion coherence for each trial [ntrials]
%       D.c: whether the choice was correct (1) or incorrect (0) [ntrials]
%
%     params: struct of parameters:
%       t: times at which to evaluate the DDM model
%       plot_flag: do (1) or do not (0) plot the data & model
%       optim_method: criteria to optimize (sum of squared
%       errors, likelihood of choice & RT, ...)
%
%
% output parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%    err: a measure of the mismatch between model and data. Depends on the
%          criterion being minimized, as specified in 'params'
%
%    P: struct with the model output.
%      P.up.p: probability of rightward choice for each unique coherence [ncoh]
%      P.up.pdf_t: probability of reaching the upper bound at each time step for 
%            each unique coherence [ncoh x ntimes]
%      P.up.cdf_t: Same as P.up.pdf_t, but cumulative over time
%      P.up.mean_t: mean decision time for each coherence for rightward choices [ncoh]
%      P.lo: same as P.up for leftward choices
%      P.drift: drift-rates for each unique coherence
%      P.Bup: upper bound
%      P.Blo: lower bound
%      P.t: vector of times


%% model parameters
kappa  = theta(1);
ndt_m  = theta(2);
ndt_s  = theta(3);
B0     = theta(4);
coh0   = theta(5);
y0a    = theta(6);
ndt_mu_delta = theta(7);

%% parse data
coh = D.coh; % motion coherence
calc_optim_criteria = true;
if isfield(D,'rt')
    rt      = D.rt; % response time
    choice  = D.choice; % choices
    c       = D.c; % correct/incorrect
else
    calc_optim_criteria = false;
end

%% make the time vector
if ~isempty(params) && isfield(params,'t')
    t = params.t;
    dt = t(2)-t(1);
else
    dt = 0.005;
    t  = 0:dt:10;
end

t = t(:);

%% plot or not?
if ~isempty(params) && isfield(params,'plot_flag')
    plot_flag = params.plot_flag;
else
    plot_flag = 0;
end

%% solve the DDM› model
Bup = B0;
drift = kappa * unique(coh + coh0);

yp = y0a/B0; % as a proportion of the bound height
P  = analytic_DDM(drift,t,Bup,yp);

%    P: struct with the model output:
%      P.up.p: probability of a rightward choice for each unique coherence [ncoh]
%      P.up.mean_t: mean decision time for each coherence for rightward choices [ncoh]
%      P.up.pdf_t: probability of reaching the upper bound at each time step for 
%            each unique coherence [ncoh x ntimes]
%      P.up.cdf_t: Same as P.up.pdf_t, but cumulative over time
%      P.lo: same as P.up for leftward choices
%      P.drift: drift-rates for each unique coherence
%      P.Bup: upper bound
%      P.Blo: lower bound
%      P.t: vector of times

%% calc. the criteria to minimize (MSE, likelihood, etc.)
if calc_optim_criteria

    if isfield(params,'optim_method')
        optim_method = params.optim_method;
    else
        optim_method = 3;
    end

    switch optim_method
        case 1 % fit choice and RT full distribution
            [err,pPred] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);


        case 2 % sum of squares errors of the mean RTs, correct trials only

            [~,mRT,eRT] = curva_media(rt,coh,c==1,0);
            udrift = P.drift;
            mean_dt = nan(size(udrift));
            mean_dt(udrift>0) = P.up.mean_t(udrift>0);
            mean_dt(udrift<0) = P.lo.mean_t(udrift<0);
            mean_dt(udrift==0) = 0.5 * (P.lo.mean_t(udrift==0) + P.up.mean_t(udrift==0)); % assumes equal 
            % number of right and left

            ndtm = nan(size(udrift));
            ndtm(udrift>0) = ndt_m + ndt_mu_delta;
            ndtm(udrift<0) = ndt_m;
            ndtm(udrift==0) = ndt_m + ndt_mu_delta/2; % assume equal number of right and left

            err = sum((mean_dt + ndtm - mRT).^2);


        case 3 % fit choice and RT full distribution, but ignore RT on incorrect trials
            [err,pPred] = logl_choiceRT_1d(P,choice,rt,coh,ndt_m,ndt_s);
            % for incorrect trials, I just take the prob. or error (ignore RT)
            ignore_RT_incorrectTrials_flag = 1;
            if (ignore_RT_incorrectTrials_flag)
                ucoh = unique(coh);
                ncoh = length(ucoh);
                pError = nan(length(ucoh),1);
                pError(ucoh>0) = P.lo.p(ucoh>0);
                pError(ucoh<0) = P.lo.p(ucoh<0);
                pError(ucoh==0) = 0.5;
                for i=1:ncoh
                    inds = c==0 & coh==ucoh(i);
                    pPred(inds) = pError(i);
                end

                %clip
                pPred(pPred<eps) = eps;

                err = -nansum(log(pPred));
            end

        case 4 % mean RT and choices, using likelihoods

            err = logl_choice_meanRT_1d(P,choice,rt,coh,c,ndt_m,ndt_m+ndt_mu_delta,ndt_s);

        case 5 % only RT means

            err = logl_choice_meanRT_1d(P,choice,rt,coh,c,ndt_m,ndt_m+ndt_mu_delta,ndt_s,'fit_choices',false);

        

    end


    %% print
    fprintf('err=%.3f kappa=%.2f ndt_mu=%.2f ndt_s=%.2f B0=%.2f coh0=%.2f y0=%.2f ndt_delta=%.2f \n',...
        err,kappa,ndt_m,ndt_s,B0,coh0,y0a,ndt_mu_delta);

else
    err = [];
end

%% plot
if plot_flag && calc_optim_criteria

    figure(1);clf
    set(gcf,'Position',[263  338  377  563])
    subplot(2,1,1);
    curva_media(choice,coh,[],3);
    hold all
    ucoh = unique(coh);
    plot(ucoh,P.up.p,'k-');
    xlim([min(ucoh),max(ucoh)]);
    xlabel('Motion coherence');
    ylabel('P rightward choice')

    subplot(2,1,2);

    % only correct trials
    rt_model_c = P.up.mean_t;
    rt_model_c(ucoh<0) = P.lo.mean_t(ucoh<0);
    rt_model_c(ucoh==0) = (P.up.mean_t(ucoh==0)+P.lo.mean_t(ucoh==0))/2;
    rt_model_c = rt_model_c + ndt_m;

    curva_media(rt,coh,c==1,3);
    hold all
    plot(ucoh,rt_model_c,'k-');
    xlim([min(ucoh),max(ucoh)]);
    xlabel('Motion coherence');
    ylabel('RT (s)')

    %set(gcf,'Position',[270   793  1084   293])

    format_figure(gcf);

    drawnow

end