function tryme
% This function identifies switching leaders using transfer entropy [3] based
% fitness function [1] in a sample dataset of 4
% particles that are simulaed using Vicsek model [2]
%

% [1] S. Butail and M. Porfiri. Detecting switching leadership in 
% collective motion. Chaos: An Interdisciplinary Journal of Nonlinear 
% Science 29, 011102, 2019
%
% [2] Vicsek, T.; Czirók, A.; Ben-Jacob, E.; Cohen, I. & Shochet, O. 
% Novel type of phase transition in a system of self-driven particles 
% Physical Review Letters, 1995, 75, 1226-1229
%
% [3] Schreiber, T. Measuring information transfer Physical Review Letters,
% 2000, 85, 461-464
% 

load('./sim_N=04_eta=3.000_etaf=3.000_ns=03_T=3200_sp=0.05_01.mat');    

% parse the signal into cells
sig=data.tr;

% data support can be based on the values obtained as in here for turn rate
% or simply [-1 1] provided you normalize the full time series as in [3]
% support=[-pi, pi];
support=[-1 1]; sig=normalize(sig);

% split into cells of size cs
cell_size=100;
if cell_size > 1
    parts=[1:cell_size:size(sig,2), size(sig,2)]; 
else
    parts=1:cell_size:size(sig,2);
end
        
for jj=1:numel(parts)-1
    if jj==numel(parts)-1
        % if this is the last part
        pp{jj}=sig(:,parts(jj):parts(jj+1));
    else
        pp{jj}=sig(:,parts(jj):parts(jj+1)-1);
    end
end

% if it can't be split into parts then just assign it as is
if numel(parts)<2
    pp{1}=sig;
end

% optimal partitioning (can be run as is to see how it works)
estm='hist'; 
bins=floor(sqrt(cell_size));
[~, pdata, ~]=opttspart(pp, ...
    @(x) gfit_teM(x, support, bins, estm));


% get the order
est.lidx=zeros(1,numel(data(1).lidx));
kk=1;
for jj=1:numel(pdata)
    npd=size(pdata{jj},2);
    [~, l0]=gfit_teM(pdata{jj}, support, bins, estm);
    est.lidx(kk+1:kk+npd)=l0;
    kk=kk+npd;
end
err=sum(est.lidx-data(1).lidx ~=0)/numel(data(1).lidx);
fprintf('err=%.3f percent \n', err*100);

figure(1); gcf; clf;
plot(data(1).lidx, 'k', 'linewidth', 4);
hold on;
plot(est.lidx, 'k--', 'linewidth', 2);

xlabel('time');
ylabel('leader id');
legend('true value', 'estimate');
set(gca, 'ylim', [1 N]);
set(gca, 'ytick', 1:N);
set(gca, 'fontsize', 16, 'fontname', 'times')
grid on
