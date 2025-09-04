function logl = logl_choice_meanRT_1d(P,choice,rt,coh,c,ndt_m0,ndt_m1,ndt_s,varargin)
% computes the logl of the parameters given the data. 
% It assumes Gaussian errors for RTs (drift-correct trials only), and binomial
% errors for the choices


fit_choices_too = true;
for i=1:length(varargin)
    if isequal(varargin{i},'fit_choices')
        fit_choices_too = varargin{i+1};
    end
end


% contribution from RT. Get the mean and standard error. Only correct
% trials.
[~,mRT,eRT] = curva_media(rt,coh,c==1,0);
udrift = P.drift;
mean_dt = nan(size(udrift));
mean_dt(udrift>0) = P.up.mean_t(udrift>0);
mean_dt(udrift<0) = P.lo.mean_t(udrift<0);
mean_dt(udrift==0) = 0.5 * (P.lo.mean_t(udrift==0) + P.up.mean_t(udrift==0)); % assume equal number of right and left

ndt_m = nan(size(udrift));
ndt_m(udrift>0) = ndt_m1;
ndt_m(udrift<0) = ndt_m0;
ndt_m(udrift==0) = (ndt_m0 + ndt_m1)/2; % assume equal number of right and left

sigma = sqrt(eRT.^2 + ndt_s.^2);
likelihood_RT = exp(-(mean_dt + ndt_m - mRT).^2./sigma.^2);

% to avoid -inf
likelihood_RT(likelihood_RT==0) = eps;

% contribution from choice
if fit_choices_too
    ucoh = nanunique(coh);

    pc = P.up.p; % pright
    likelihood_Choice = nan(size(ucoh));
    for i=1:length(ucoh)
        I = coh==ucoh(i);
        nc = sum(choice(I)==1);
        likelihood_Choice(i) = binopdf(nc,sum(I),pc(i));
    end
    % to avoid -inf
    likelihood_Choice(likelihood_Choice==0) = eps;
    
    logl = -1 * (nansum(log(likelihood_RT)) + nansum(log(likelihood_Choice)));
    
else
    
    logl = -1 * nansum(log(likelihood_RT));
    
end





