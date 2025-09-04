function [choice, DecT, DV] = bound_cross_flat(times, dec_var, bound)
% Detects the crossing of a flat-bound in simulated data. 

% INPUTS: 
% times: time vector [num_times]
% dec_var: decision variable [num_trials x num_times]
% bound: upper bound height [scalar value]

% OUTPUTS: 
% choice: 0 or 1 if the lower or upper bound is reached. It is NaN is no
% bound is crossed [num_trials]
% DecT: time taken by the decision variable to reach a bound. It is NaN if
% no bound was reached[num_trials]
% DV: same as dec_var, except that the values after the crossing of the
% bound are filled with NaNs. Optional output. [num_trials x num_times]. 


if nargout>2
    output_trunc_dv = 1;
else
    output_trunc_dv = 0;
end


[ntrials, ntimes] = size(dec_var);


% add a last column with forced crossing
dv_up = [dec_var, ones(ntrials,1)*bound*1.1];
dec_step_up = single_bound_cross_dynamic(dv_up, bound);

dv_lo = [-1*dec_var, ones(ntrials,1)*bound*1.1];
dec_step_lo = single_bound_cross_dynamic(dv_lo, bound);

[y,isWinner,val] = winning_race([dec_step_up, dec_step_lo]);

no_crossing = sum(isWinner,2)>1; % identify the trials with no crossing
choice = double(y==1);
choice(no_crossing)=nan;

dec_step = val;
DecT = nan(ntrials,1);
I = dec_step<=ntimes;
DecT(I) = times(dec_step(I));


if output_trunc_dv
    DV = dec_var;
    for i=1:ntrials
        DV(i,(dec_step(i)+1):end) = nan;
    end
end


end


function t_step = single_bound_cross_dynamic(dv, Bup)

a = dv>Bup';
t_step = sum(cumprod(~a,2),2)+1;
nt = size(dv,2);
I = t_step==(nt+1);
t_step(I) = nan;

end

function [y,isWinner,val] = winning_race(dt)
% function y = winning_race(dt)
%dt [ntrials x races]: decision time of each race. NaN means
%bound was not reached
%y [ntrials x 1]: id of winner. If two are equal, returns the first
%isWinner: to determine ties.

[val,y]       = min(dt,[],2);
y(isnan(val)) = nan;
isWinner      = bsxfun(@eq,dt,val);
end