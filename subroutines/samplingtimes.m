function [Times] = samplingtimes(ZZ,no_samples,which_gridfunction)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%rng('shuffle')% add as optional input default 'shuffle'
largeT=1e6/std(ZZ(:).^2);% add as optional input
smallT=1e-12*largeT; % add as optional input

[no_snps,no_traits]=size(ZZ);
rows=randsample(1:no_snps,no_samples,true,sum(ZZ.^2,2).*prod(ZZ~=0,2));
Times=1./ZZ(rows,:);

%
Times(Times>largeT)=largeT;Times(Times<-largeT)=-largeT; % prevents overflow issues

% Add in times close to the origin
if ~exist('which_gridfunction')
    which_gridfunction=1;
end
if which_gridfunction==1
    temp = get_smallgrid(smallT,no_traits);
elseif which_gridfunction==2
    temp = get_grid([-smallT 0 smallT]',no_traits);
else
    temp=[];
end
Times=[Times;temp];

    function [grid] = get_grid(gridvals,no_traits)
        
        
        gridsize=length(gridvals);
        no_combos=gridsize^(no_traits);
        grid=zeros(no_combos,no_traits);
        for kk=1:no_traits
            row=repmat(gridvals,gridsize^(no_traits-kk),gridsize^(kk-1))';
            grid(:,kk)=row(:);
        end
    end

    function [grid] = get_smallgrid(smallT,no_traits)
        
        
        pairs=nchoosek(1:no_traits,2);
        grid=zeros(nchoosek(no_traits,2)*4+no_traits*2+1,no_traits);
        for ii=1:size(pairs,1)
            grid(ii*4-3:ii*4,pairs(ii,:))=smallT*[1 1; 1 -1; -1 1; -1 -1]/sqrt(2);
        end
        for ii=1:no_traits
            grid(size(pairs,1)*4+(ii*2-1:ii*2),ii)=smallT*[1; -1];
        end
        
    end


end

