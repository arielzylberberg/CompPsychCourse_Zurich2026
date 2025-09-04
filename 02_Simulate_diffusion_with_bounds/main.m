addpath(genpath('../functions_addtopath/'));

%%

% seed = 104046;
% rng(seed,'twister');

%% some figure defaults
set(0,'defaultAxesFontSize',18);
set(0, 'DefaultLineLineWidth', 1);

%%  parameters

% parameters
NumTrials    = 5000; % number of simulated trials
MaxTime      = 5; 
dt           = 0.0001; % time step
sigma        = 1; % standard deviation of the DV after 1 sec. of accumulation
mu           = 1; % mean of the DV after 1 sec of accumulation
T = 0:dt:MaxTime; % time vector

% momentary evidence
e = mu * dt + sigma * sqrt(dt) * randn(NumTrials, length(T));

% accumulated evidence
dv = cumsum(e,2);


%%

Bup = 1;
[choice, DecT, DVnan] = bound_cross_flat(T, dv, Bup);

%% plots
% some trajectories
p = publish_plot(1,1);
plot(T, DVnan(1:100,:)');
xlabel('Time [s]');
ylabel('Accumulated evidence');
p.format();

% distribution of the DV for a fixed time
p = publish_plot(1,1);
histogram(DVnan(:,5000),100);
xlabel('Decision variable (a.u.)');
ylabel('# trials');
p.format();

% decision time vs. accuracy
p = publish_plot(1,1);
[tt,xx,ss] = curva_media(DecT, choice, [],0);
terrorbar(tt,xx,ss);
xlabel('Accuracy');
ylabel('Decime time [s]');
xlim([-0.1,1.1]);
p.format();

p = publish_plot(1,1);
histogram(DecT,100);
ylabel('# trials');
xlabel('Decision Time [s]');
p.format();


