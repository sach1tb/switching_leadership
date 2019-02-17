function h=ent(data, numberOfBins, support, form, correction)
% function h=ent(data, numberOfBins, support, form, correction)

% ent(X, numberOfBins, support) assumes that you have multiple values of a variable and you
% want to bin it into numberOfBins bins. p doesn't mean p but a variable

h=0;
X=data(1,:); 
if size(data,1)>1, Y=data(2,:); end

switch form
    case 'x' % a single variable entropy
        bins=linspace(support(1), support(2), numberOfBins);

        % histograms
        p=hist(X, bins);
        p=p/sum(p);
%         p(p==0)=1;
        p=p(p~=0);
        h=-sum(p(:).*log2(p(:)));
        
        % Miller-Madow correction
        if nargin >4
            if strcmp(correction, 'millermadow')
                h=h+(numel(p)-1)/(2*size(X,2));
            end
        end
        
    case 'x;y' % joint entropy
        for jj=1:2
            bins{jj}=linspace(support(jj,1), support(jj,2), numberOfBins(jj));
            binwidth=bins{jj}(2)-bins{jj}(1);
            bb=bins{jj};
            bins{jj}=[-inf bins{jj}(1)-binwidth/2 bins{jj}+binwidth/2];
            bins{jj}(end)=bb(end);
            nbins(jj)=numel(bins{jj});
        end
        % p(X[markovOrder-1], Y[markovOrder-1])
        pXY=zeros(nbins(1)-1, nbins(2)-1);
        for jj=1:nbins(1)-1
            for kk=1:nbins(2)-1
                pXY(jj,kk)=sum(X>=bins{1}(jj) & X<=bins{1}(jj+1) & ...
                                       Y>=bins{2}(kk) & Y<=bins{2}(kk+1));
            end
        end
        pXY=pXY(2:end,2:end);
        
        pXY=pXY(pXY~=0);

        
        pXY=pXY/sum(pXY(:));
                
        h=-sum(pXY(:).*log2(pXY(:)));
        
    case 'x|y' % conditional entropy
        
        bins=linspace(support(1), support(2), numberOfBins);
        binwidth=bins(2)-bins(1);
        bb=bins;
        bins=[-inf bins(1)-binwidth/2 bins+binwidth/2];
        bins(end)=bb(end);
        nbins=numel(bins);
        
        % p(X[markovOrder-1], Y[markovOrder-1])
        pXY=zeros(nbins-1, nbins-1);
        for jj=1:nbins-1
            for kk=1:nbins-1
                pXY(jj,kk)=sum(X>=bins(jj) & X<=bins(jj+1) & ...
                                       Y>=bins(kk) & Y<=bins(kk+1));
            end
        end
        pXY=pXY(2:end,2:end);
        pXY=pXY/sum(pXY(:));
        
        % pY
        pY=zeros(nbins-1,1);
        for jj=1:nbins-1
            pY(jj)=sum(Y>=bins(jj) & Y<=bins(jj+1));
        end
        pY=pY(2:end);
        pY=pY/sum(pY(:));
        
        for ii=1:nbins-2
            for jj=1:nbins-2
                if pXY(ii,jj) && pY(ii)
                    h=h-pXY(ii,jj)*log2(pXY(ii,jj)/pY(ii));
                end
            end
        end    
end

