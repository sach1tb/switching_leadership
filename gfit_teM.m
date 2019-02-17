function [tmax, idx]=gfit_teM(data, support, bins, estm)
% function [tmax, idx]=gfit_teM(data, support, bins, estm)
%
% transfer entropy fitness function
%
%

ns=size(data,2);
N=size(data,1);

% adaptive binning
% if ~bins(1)
%     bins(1)=min(floor(sqrt(ns)), 10);
% end
TEi2j=zeros(N,N,numel(bins));

for bb=1:numel(bins)
    for ii=1:N
        for jj=1:N
            if ii~=jj
               TEi2j(ii,jj,bb)=ete(estm, data(jj,:), data(ii,:), 1, ...
                   bins(bb), 1/bins(bb), support); 
               
               % normalize w.r.t. square root of entropy of each
               Hjj_kp1=ent(data(jj,2:end), bins(bb), support, 'x');
               Hii_k=ent(data(ii,1:end-1), bins(bb), support, 'x');
               if Hjj_kp1*Hii_k > eps
                   TEi2j(ii,jj,bb)=TEi2j(ii,jj,bb)/sqrt(Hjj_kp1*Hii_k);
               else
                   TEi2j(ii,jj,bb)=0;
               end
            end
        end
    end
end

% find net TE
netTEi2j=TEi2j-TEi2j';

[tmax, idx]=max(sum(sum(netTEi2j,3),2)/(N-1)/(ns^(-1/2)));