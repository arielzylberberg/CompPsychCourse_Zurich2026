addpath('../functions_addtopath/');

% parameters
NumTrials    = 5000;
MaxTime      = 1;
dt           = 0.0001;
sigma        = 1;
mu           = 1;
T = 0:dt:MaxTime;

% momentary evidence
e = mu * dt + sigma * sqrt(dt) * randn(NumTrials, length(T));

% accumulated evidence
dv = cumsum(e,2);

% plot 
p = publish_plot(1,1);
plot(T, dv(1:100,:)');
xlabel('Time [seconds]');
ylabel('Decision variable [a.u.]');
p.format('Presentation');

% mean and variance of accumulated evidence
p = publish_plot(2,1);
p.next();
plot(T, mean(dv));
ylabel('Mean of the DV');

p.next();
plot(T, var(dv));
ylabel('Variance of the DV');
xlabel('Time [s]');
p.format('Presentation');

