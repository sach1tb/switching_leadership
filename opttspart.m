function [pmax, pdata, optval]=opttspart(data, g)
% opttspart optimally partitions time series data into successive blocks 
% that maximize a fitness function g
%
% Syntax
%
%   opttspart
%   [pmax, pdata, optval]=opttspart(data, g)
%
% Description
%   opttspart as is will run example test cases on a square wave with varying 
%   frequency. An additive function g=@(x) max(sum(-x), sum(x)) picks out 
%   the points where the values change from positive to negative and vice versa.
%   Two test cases can be run by setting the variable tc=1 or 2 
%   tc=1 will partition the square wave without cells tc=2 will partition the
%   same square wave but instead with the time series sectioned into
%   equidistant cells of size 10.
%
%   [pmax, pdata, optval]=opttspart(data, g) takes as input an m x n data 
%   matrix of m time series each of length n steps (or cells) and a fitness 
%   function g for optimally partitioning [1]. The output pmax is the location 
%   of the partition (time step or cell number), pdata is the time series 
%   data partitioned into cells according to pmax and optval is the optimal 
%   value g*; g(0) shown on the title in the example is the value of the 
%   fitness function with no partitions
%
%   %% Using cells [2] 
%   The algorithm has O(n^2) operations and therefore depending on the fitness
%   function, can be slow. In [2] we use cells section the time series into 
%   cells before optimal partitioning. This makes the algorithm faster but
%   introduces an error if a change does not occur exactly where two cells are
%   sectioned.
% 
% References
%   [1] Jackson, B., Scargle, J. D., Barnes, D., Arabhi, S., Alt, A., Gioumousis,
%   P., ... & Tsai, T. T. (2005). An algorithm for optimal partitioning of 
%   data on an interval. IEEE Signal Processing Letters, 12(2), 105-108. This
%   is the paper with the algorithm coded here. 
%
%   [2] Butail, S. and Porfiri, M. Detecting switching leadership in collective
%   motion, Chaos, 29, 011102, 2019. This paper uses this algorithm with an
%   information theoretic fitness function to infer causality in pairs of
%   time series. 
%
%


%
% ** general algorithm flow 
%
% 1) opt(0)=0
% 2) j=0,0
% % first for loop
% 	end(1,1)=g(B1:B1) % j is 1 because there are two conditions j goes 
%   from 1-n+1 and n is the maximum value available for j, so end(j, n+1) 
%   is end (1, 0+1)
%
% 	[opt(1), idx]=max(opt(0), end(1,1))
% 	lastchange(1)=idx
% 
% 	end(2,2)=g(B1:B2) % second iteration and 1+1
% 	[opt(2), idx]=max(opt(0)+end(1,1), opt(1)+end(2,2))
% 	lastchange(2)=idx
% 
% 	end(3,3)=g(B1:B3)
% 	[opt(3), idx]=max(opt(0)+end(1,1), opt(1)+end(2,2), opt(2)+end(3,3))
% 	lastchange(3)=idx
% 
% ...
% 	end(T,T)=g(B1:BN)
% 	[opt(T), idx]=max(...)
% 	lastchange(T)=idx

% test cases
tc=1;

if nargin < 1
    t=0:0.1:12*pi;
    data=sign(sin(2*t.^.6));
    if tc==2, cs=10; else cs=1; end

    % maximum value between the sum of elements or -ve sum of elements if
    % the partition has only -ve elements then a -ve of that is more and if
    % it has a mix of -ve and +ve then it is better to move it until you
    % get a -ve or +ve only
    g=@(x) max(sum(-x), sum(x));

    if cs > 1
        parts=[1:cs:size(data,2), size(data,2)];
    else
        parts=1:cs:size(data,2);
    end
    for jj=1:numel(parts)-1
        dataf{jj}=data(:,parts(jj):parts(jj+1)-1);
    end
    data=dataf;
end

T=size(data,2);
opt(1)=0;

% O(T^2) operations to compute end1
% define end(j,n) = g(Bj:n)
% we now have a matrix where each entry corresponds to the fitness level of
% the partition from any point to the end point n
for n=1:T
    for j=1:n
        end1(j,n)=g(cat1(data, j,n));
    end
end

% this loop creates a new array of all possible partitions as below
% [ opt(1) + g(B2:n)
%   opt(2) + g(B3:n)
%   opt(3) + g(B4:n)
%   ...
%   opt(n-1) + g(Bn-1:n)]
% and then selects the best partition in terms of fitness
% it then stores that index in a lastchange
% 
lastchange(1)=0;
for n=2:T
    arr=[];
    for j=2:n
        arr=[arr; opt(j-1)+end1(j,n)];
    end
    [opt(j), lastchange(j)]=max(arr);
end

% this loops through lastchange to select the partitions
ni=T;
pmax=lastchange(T);
while ni >1
    ni=lastchange(ni-1);
    pmax=[pmax, ni];
end
% remove the first zero
pmax=pmax(pmax~=0);

% also throw the first partition if there are more
% instead make sure to say that pmax starts at 1
if numel(unique(pmax)) > 1
    pmax=pmax(pmax~=1);
end

pmax=unique(sort(pmax));

% divide the data into partitions
np=numel(pmax);
if np >=2
pdata{1}=cat1(data, 1, pmax(1));
for ii=2:np
    pdata{ii}=cat1(data, pmax(ii-1)+1, pmax(ii));
end
pdata{ii+1}=cat1(data, pmax(ii)+1, T);
else
    pdata{1}=cat1(data, 1, T);
end
optval=opt(j);


% show figure for example cases; no figure is plotted otherwise
if nargin < 1 
    if iscell(data)
        data=cat(2, data{:});
    end
    figure(1); gcf; clf;
    plot(1:numel(data), data, 'k', 'linewidth', 2);
    hold on;
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(ones(2,1)*pmax*cs, [min(data); max(data)]*ones(1,numel(pmax)), ...
        'r-', 'linewidth', 2);
    xlabel('time');
    ylabel('data');
    title(sprintf('g^*=%.3f v/s g(0) = %.3f', opt(j), ...
        g(data(1:T))));
    grid on
    legend('time series', 'partitions');
end

% common function to concatenate cells and values 
function data=cat1(data, n1, n2)
if iscell(data)
    data=cat(2, data{n1:n2});
else
    data=data(:,n1:n2);
end
