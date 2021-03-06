function [phiB, P, noiseVarRemaining, times] = ecfMake(data,times,varargin)
% ecfMake creates the fields needed for the ECF object by calculating
%   sampling times, projection matrix, phi values, and projection matrix to
%   prune highly correlated sampling times.
%
% Required Inputs:
%   data: DATA object containing GWAS summary statistics
%   times: Vector of sampling times
%
% Optional Inputs:
%   whichTraits: Which traits from data to include. Can be vector of indices,
%       or cell array of trait names. Order matters.
%   maxNoSamplingtimes: Maximum number of sampling times
%   noBlocks: Number of jackknife blocks
%   tol: % Tolerance for projection matrix. If this value is too small,
%       seems to produce miscalibrated p-values
%
% Outputs:
%   phiB: phi values for each jackknife block and sampling time combination
%   P: Projection matrix for pruning highly correlated sampling time
%       combinations
%   noiseVarRemaining: Covariance matrix of the noise in the data
%   times: (# traits x # sampling time combinations) array of all sampling
%       time combinations

p=inputParser;

addRequired(p, 'data', @(obj)isa(obj,'DATA'));
addRequired(p, 'times', @ismatrix);

addParameter(p, 'whichTraits', []);
addParameter(p, 'maxNoSamplingtimes', 1e3);
addParameter(p, 'noBlocks', 100);
addParameter(p, 'tol', 0);

parse(p,data,times,varargin{:});

% Identify which traits to use
whichTraits=p.Results.whichTraits;
if isempty(whichTraits)
    traitIdx=1:data.noTraits;
end

if iscell(whichTraits)
    [~,i1,traitIdx]=intersect(whichTraits,data.traits,'stable');
    if length(i1)~=length(whichTraits)
        error('Specified traits not found in data')
    end
end

maxNoSamplingtimes = p.Results.maxNoSamplingtimes;
times=get_grid(times',length(traitIdx));
if size(times,2) > maxNoSamplingtimes
    incl = randsample(1:size(times,2),maxNoSamplingtimes,false);
    times=times(:,incl);
end

% Construct jackknife blocks
if isempty(data.blocks)
    noBlocks=p.Results.noBlocks;
    blocksize=floor(data.noSNPs/noBlocks);
    data.blocks=arrayfun(@(b){(b-1)*blocksize+1:min(b*blocksize,data.noSNPs)},(1:noBlocks)');
else
    noBlocks=length(data.blocks);
end

noiseVar=data.sigmaEps(traitIdx,traitIdx);

Z = data.z;

ecf_fn=@(z,noiseScalars)mean(cos(z*times),1)-1;

% ECF of each jackknife block
phiB=cellfun(@(i){ecf_fn(Z(i,traitIdx),1)},data.blocks);

phiB=vertcat(phiB{:})';

if p.Results.tol > 0
    % normalize rows of phiB
    D=std(phiB,[],2);
    
    % Truncated-svd regularization
    [U,S,~]=svd((phiB-mean(phiB,2))./D,'econ');
    
    % Projection matrix whitens residuals + truncates small eigenvalues
    S=diag(S);
    incl=S > mean(S) * p.Results.tol;
    Sinv=diag(1./S(incl));
    
    P=Sinv*U(:,incl)'./D';
    
    fprintf('ECF contains %d combinations of sampling times\n',sum(incl))
    if sum(incl)<10
        warning('Very few combinations of sampling times are retained; use more jackknife blocks or a smaller tolerance parameter')
    end
else
    P = eye(size(times,2));
end

% noiseVar if not already corrected out
noiseVarRemaining = noiseVar;
end