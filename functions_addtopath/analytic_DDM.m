function D =  analytic_DDM(drift,t,Bup,yp,ny)
% D =  analytic_DDM(drift,t,Bup,yp,ny)
% Implementation of Cox & Miller (PDF of non-absorbed) and Ratcliff
% (probability absorbed and first passage time) analytic solutions to bounded
% drift diffusion (drift to bound). The solutions here handle flat but
% asymmetric bounds with delta function initiation anywhere between bounds.
%
% Input arguments
% ~~~~~~~~~~~~~~~
% drift         vector of drift rates
% t             time series or set of times in seconds (t(1) needs to be zero)
% Bup           scalar bounds with   Blo =-Bup (for asymmertric bounds
%               change yp instead
% yp            value of DV at t=0 (scalar) as a poprtion of Bup [i.e, -1 to +1]
% ny (optional) specifies the granularity (number of points) of  DV from [Blo,Bup]
%               and asks the function to return the unabsorbed probability distribution 
%               all times (will slow down simulation).
%
% Outputs
% ~~~~~~~~~~~~~~~
% Returns D, a structure (the first four have "lo" vesions too)
% D.up.p(drift)        probability of hitting the upper bound for each drift
% D.up.mean_t(drift)   mean decision time for upper bound for each drift level
% D.up.pdf_t(t,drift)  probability of upper bound hit at each time (sums to D.up)
% D.up.cdf_t(t,drift)  cumulative probability of upper bound hit at each time (ends at D.up)
%
% D.drifts             drifts tested
% D.bounds             bounds
% D.t                  times used for simulation
%
% if ny is passed then th:e function also returns:
% D.notabs.pdf(drifts,y,t) probability of being at y at time t
% D.notabs.pos_t(drifts,t) probability of not being absorbed and being >0
% D.notabs.neg_t(drifts,t) probability of not being absorbed and being <0
% D.notabs.y               the sets of y's considered for the pdf
%
% This is a beta version from Daniel Wolpert.
%
% History 10/2011 Daniel Wolpert wrote it
%         26/10/2011 mns added comments
% v1.01 Jan 2nd 2013 fixed the zero coherence condition which was wrong for non symmetric
% bounds
% 07/2019 Parallelized by Ariel Zylberberg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Ratcliff formulation lower bound needs to be at 0, so set upper at B and adjust y0
yp = clip(yp,-1,1);

Blo =  -Bup;
y0  = Bup*yp;
B   = Bup-Blo;
y0  = y0-Blo;


%to calculate lower crossings information [Ratcliff actually calculates lower]
if nargin ==5
    [D.lo D.notabs]=run_analytic(t,+drift,y0,B,ny);
else
    [D.lo D.notabs]=run_analytic(t,+drift,y0,B);
end

%to calculate upper crossings (flip drift and adjust starting point)
D.up = run_analytic(t,-drift,B-y0,B);

D.t = t;
D.drift = drift;
D.Bup = Bup;
D.Blo = Blo;

%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D E] = run_analytic(t,drift,y0,B,ny)

% this uses A8 and A12 equation for Ratcliff psych review 1978

ndrift=length(drift);
nt=length(t);

% test if the user puts in a time series or an arbitrary set of times
if nt>1 && var(diff(t))<100*eps
    series_flag=1;
else
    series_flag=0;
end

%G is survivor function for zero bound
G=zeros(ndrift,nt);

%analytic expression for probability of crossing upper bound from Ratcliff A8
P=(exp(-2*drift*B)-exp(-2*drift*y0))./(exp(-2*drift*B)-1);

P(abs(drift)<100*eps)=1-y0/B; %fix for zero drift condition.

%this sums over 1->Inf but we go 1:K(end), this seems to sum far enough for standard parameters
K=1:150;

p_threshold=1e-8; %threshold for proportion un-terminated - used for series only

% this next part uses A12 equation from Ratcliff psych review 1978

for j=1:ndrift
    d=(pi/B^2)*exp(-y0*drift(j));
    den=(drift(j)^2 + pi^2*K.^2/B^2);

    
    aa = exp(-0.5*(drift(j)^2+pi^2*K.^2/B^2)'*t');
    cc = 2*K.*sin(K*pi*y0/B);
    NUM = bsxfun(@times,cc',aa);
    
    ss = bsxfun(@times,NUM',1./den);
    ss = sum(ss,2);
    G(j,:) = d*ss;

%     for k=1:nt
%         num=2*K.*sin(K*pi*y0/B).*exp(-0.5*(drift(j)^2+pi^2*K.^2/B^2)*t(k));
% 
%         ss=sum(num./den);
%         
%         G(j,k)=d*ss;   % proportion left still to cross upper bound
%         
%         %stop if we are doing a time series and are below threshold for
%         %un-terminated after at least 10 steps
%         if series_flag && k>10 &&  G(j,k)<p_threshold ,break, end
%     end
    
    
    %     % provide warning to increase maximum time if not below threshold  by maximum time
    %     if series_flag && G(j,k)>p_threshold && 0
    %         fprintf('Warning: at max time [%.2fs] %.5f%% yet to hit bound for drift=%f\n',...
    %             max(t),G(j,k)*100,drift(j));
    %     end
    G(j,:) = P(j)-G(j,:);
end

%series does not converge properly for t=0 but we know that nothing has crossed
%bound so set to 0 for all drifts if t==0
G(:,t==0)=0;


if series_flag
    %transform cumulative into pdf of stopping times and pad beginning with zeros
    dtdist=[zeros(ndrift,1) diff(G')' ]'; % this assumes t(0)=0;
    
    %average termination times (subtract half sampling interval as we used diff)
    dts=t'*dtdist;
    dts=dts./P'-(t(2)-t(1))/2;
else
    dts=[];
    dtdist=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%from cox and miller eqn 78 page 222 to get full pdf of DV between bounds

if nargin ==5 %only calculate if we received ny
    a=B-y0;
    b=y0;
   
    y=linspace(-b,a,ny)';
    [Y,T]=meshgrid(y,t+eps); % speed up by using matrices for evaluations
    a1=1./(sqrt(2*pi*T));

    %the series_flag sums over -Inf->Inf but we go -nk->nk and this seems far enough for standard parameters
    nk=10;
    
    for j=1:ndrift
        pdf=zeros(nt,ny);
        
        for i=-nk:nk
            yn=2*i*(a+b);
            ynn=2*a-yn;
            
            a2=exp(drift(j)* yn - (Y- yn-drift(j)*T).^2./(2*T));
            a3=exp(drift(j)*ynn - (Y-ynn-drift(j)*T).^2./(2*T));
            
            pdf= pdf+a1.*(a2-a3);
        end
        
        pdf=pdf*(y(2)-y(1));
        
        if t(1)==0
            pdf(1,:)=0;
            pdf(1,round(ny/2))=1;
        end
        E.pdf(j,:,:)=pdf;
        E.pos_t(:,j)=sum(pdf(:,y>=0),2);
        E.neg_t(:,j)=sum(pdf(:,y<0),2);
    end
    E.y=y;
else
    E.pdf=[];
    E.y=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D.p=P;
% D.cdf_t=G';
D.cdf_t = G;
%D.mean_t=dts;
D.mean_t=(B.*coth(B.*drift) - coth(drift.*(B - y0)).*(B - y0))./drift;
D.mean_t(drift==0)=(B/2)^2;

% D.pdf_t=dtdist;
D.pdf_t=dtdist';


