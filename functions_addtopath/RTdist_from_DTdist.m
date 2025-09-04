function rt_pdf = RTdist_from_DTdist(t, dect_pdf, ndt_m, ndt_s)
% function RTdist_from_DTdist(dect_pdf, ndt_m, ndt_s)

% Convolves the decision time distribution with the distribution of
% non-decision times

pd = makedist('Normal','mu',ndt_m,'sigma',ndt_s);
pd_trunc = truncate(pd,0,inf);
dt = t(2)-t(1);
ndt = pd_trunc.pdf(t)*dt;

% should speed things up
imax = find(cumsum(ndt)>0.999999,1);
ndt = ndt(1:imax);

rt_pdf = conv2(1,ndt(:),dect_pdf);
nt = length(t);
rt_pdf = rt_pdf(1:nt);


