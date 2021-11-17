function theta = sample_datacov(datacov,cpt_type,varargin)
% Generates theta for either a full-rank or rank-one component by sampling
% from multivariate normal distribution with mean 0 and covariance = datacov

p=inputParser;

addRequired(p, 'datacov');
addRequired(p, 'cpt_type');

addParameter(p, 'nsamples', size(datacov,1));

parse(p,datacov,cpt_type,varargin{:});

nsamples=p.Results.nsamples;


no_traits = size(datacov,1);
sigma = zeros(no_traits);

if cpt_type == "full_rank"
    for p=1:nsamples
        x = mvnrnd(zeros(no_traits,1),datacov,1);
        sigma = sigma + x'*x;
    end
    
    sigma = sigma./nsamples;
    if det(sigma) < 1e-15
        sigma = sigma + 1e-12*eye(length(sigma));
    end
    theta=triuind(chol(sigma));
    
elseif cpt_type == "rank_one"
    theta = mvnrnd(zeros(no_traits,1),datacov,1)';
    
else
    error("cpt_type must be either 'full_rank' or 'rank_one'");
    
end

end

