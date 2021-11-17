function [f, df, theta, thetaCon] = component_function(name, noTraits, fixedParam, theta)
%CptFn returns a function handle for the specified type of component
%   name: componet type. noTraits: number of phenotyeps. whichTrait: eg for
%   the trait_specific component type, index of that trait.

cpt_names={'trait_specific','rank_one','null_cpt','uncorr','uniform_corr','uniform_corr_slope','full_rank','fixed','trait_subset'};

[~,~,i]=intersect(name,cpt_names);

if isempty(i)
    error('Component name invalid')
end

% pick random theta
if nargout > 2 || nargin < 4
    randomThetaFn={@()randn(1)^2;
        @()randn(noTraits,1);
        @()[];
        @()randn(noTraits,1).^2;
        @()rand(1)*2-1;
        @()[rand(1)*2-1; randn(noTraits-1,1).^2];
        @()randn(noTraits*(noTraits+1)/2,1);
        @()[];
        @()randn(length(fixedParam)*(length(fixedParam)+1)/2,1)};
        
    theta=randomThetaFn{i}();
    
    thetaCon={[0 inf];
        repmat([-inf inf],noTraits,1);
        zeros(0,2);
        repmat([0 inf],noTraits,1);
        [-1 1];
        [-1 1; repmat([0 inf],noTraits-1,1)];
        repmat([-inf,inf],noTraits*(noTraits+1)/2,1);
        zeros(0,2);
        zeros(0,2)
        };
    thetaCon=thetaCon{i};
    
 
end

[f,df]=feval(name,theta,noTraits,fixedParam);




end

