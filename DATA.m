classdef DATA < matlab.mixin.Copyable
    properties
        traits {iscell} % trait names or file names
        z {ismatrix} % Z score or effect-size estimate matrix
        beta {ismatrix} % true causal effect sizes for simulated data
        alpha {ismatrix} % true marginal effect sizes for simulated data
        snps {isvector} % vector of truncated RSIDs or other identifier
        blocks {iscell} % LD blocks
        sigmaEps {ismatrix} % covariance matrix of sampling noise
        LDeigenvectors % cell array of eigenvectors for LD blocks
        LDeigenvalues % cell array of eigenvalues
        l2
        weights
        ldsc
        additive_cpts
        whichSNPs % cellarray of lists of which SNPs belong to each component
        normalizer
        noSNPsAnnot % number of SNPs in each l2 annotation
    end
    
    properties (Dependent)
        noTraits
        noSNPs
        noBlocks
        cov
    end
    
    methods
        % class constructor
        function obj = DATA(input,varargin)
            
            if nargin == 0
                return;
            end
            
            p=inputParser;
            % Primary input
            addRequired(p, 'input', @(obj)isa(obj,'MODEL') || iscell(obj) || ischar(obj));
            
            % Number of SNPs
            addParameter(p, 'mm', []);
            
            % Number of jackknife blocks
            addParameter(p, 'noBlocks', 100);
            
            % Eigenvectors of LD matrices (one cell per block)
            addParameter(p, 'LDeigenvectors', {});
            
            % Eigenvalues of LD matrices (one column per block)
            addParameter(p, 'LDeigenvalues', []);
            
            % Path to LD scores
            addParameter(p, 'LDscoreDir', []);
            % Path to regression weights
            addParameter(p, 'WeightsDir', []);
            
            % pre-loaded LD scores
            addParameter(p, 'l2snp', []);
            addParameter(p, 'l2', []);
            
            % Covariance matrix of epsilon
            addParameter(p, 'noiseVar', []);
            
            % Whether to use additive components
            addParameter(p, 'additive_cpts', 1);
            
            parse(p,input,varargin{:});
            
            % load LD data
            if ~isempty(p.Results.LDeigenvalues)
                noBlocks = numel(p.Results.LDeigenvectors);
                if any(size(p.Results.LDeigenvalues) ~= size(p.Results.LDeigenvectors))
                    error('LDeigenvalues and LDeigenvectors should be cell arrays of the same size')
                end
                
                % discard eigenvalues equal to zero or NaN (+ associated
                % eigenvectors)
                obj.LDeigenvalues = p.Results.LDeigenvalues;
                obj.LDeigenvectors = cellfun(@(A,S){A(:,1:length(S))},...
                    p.Results.LDeigenvectors, obj.LDeigenvalues);
                blocks = [0,cumsum(cellfun(@length,obj.LDeigenvectors))];
                mm=blocks(end);
                obj.blocks = arrayfun(@(first,last){first+1:last},blocks(1:end-1),blocks(2:end));
            else
                mm=p.Results.mm;
            end
            
            obj.additive_cpts = p.Results.additive_cpts;
            
            noisevar = p.Results.noiseVar;
                
            
            % simulate from MODEL
            if isa(input,'MODEL')
                if isempty(mm); error('Please specify number of SNPs to simulate from a MODEL'); end
                
                [obj.beta, obj.whichSNPs] = input.simulate(mm);
                obj.sigmaEps = noisevar;
                
                % simulate sumstats with LD
                if ~isempty(obj.LDeigenvalues)
                    % R*beta
                    alpha = cellfun(@(i,U,S){U*diag(S)*U'*obj.beta(i,:)}...
                        ,obj.blocks, obj.LDeigenvectors, obj.LDeigenvalues);
                    
                    % noise ~ MVN with 2 covariance matrices. Row
                    % covariance (across SNPs) is proportional to R; column
                    % covariance (across traits) is equal to noisevar
                    noise = cellfun(@(U,S){U*diag(sqrt(S))*...
                        mvnrnd(zeros(1,input.noTraits),noisevar,length(S))}...
                        ,obj.LDeigenvectors,obj.LDeigenvalues);
                    
                    obj.alpha = vertcat(alpha{:});
                    obj.z = vertcat(alpha{:}) + vertcat(noise{:});
                    
                    % simulate sumstats without LD
                else
                    obj.alpha = obj.beta;
                    obj.z = obj.beta + mvnrnd(zeros(1,input.noTraits),noisevar,mm);
                end
                
                obj.traits = arrayfun(@(i){['trait ',num2str(i)]},1:input.noTraits);
                
                % import files on paths specified in cell array
            elseif (iscell(input) || ischar(input))
                if iscell(input)
                    [obj.snps,obj.z,phase]=import_sumstat_files(input);
                    
                    % truncate .sumstats from input
                    obj.traits=cellfun(@(x)x(1:end-9),input,'UniformOutput',false);
                    
                    % import files in directory specified as string
                elseif ischar(input)
                    % find files in directory
                    files=dir([input,'/*.sumstats']);
                    if isempty(files)
                        error('Input string should point to a directory containing sumstats files with .sumstats extension')
                    end
                    
                    % import data from those files
                    [obj.snps,obj.z]=import_sumstat_files(...
                        cellfun(@(s){[input,'/',s]},{files.name}));
                    
                    % truncate .sumstats from filenames
                    obj.traits=cellfun(@(x){x(1:end-9)},{files.name});
                end
                
                l2snp = p.Results.l2snp;
                l2 = p.Results.l2;
                if ~isempty(p.Results.LDscoreDir) || (~isempty(l2snp) && ~isempty(l2))
                    disp('Importing LD scores')
                    if isempty(l2snp) || isempty(l2)
                        [l2snp,l2,column_names,noSNPsAnnot] = load_l2file(p.Results.LDscoreDir,'.l2.ldscore');
                    end
                    [wtsnp,weights] = load_l2file(p.Results.WeightsDir,'.l2.ldscore');
                    [l2snp,l2_idx,wt_idx] = intersect(l2snp,wtsnp,'stable');
                    l2=l2(l2_idx,:);
                    weights=weights(wt_idx);
                    obj.noSNPsAnnot = noSNPsAnnot;
                    
                    [obj.snps, l2_idx, z_idx] = intersect(l2snp,obj.snps,'stable');
                    obj.z = obj.z(z_idx,:);
                    obj.l2 = l2(l2_idx,:);
                    obj.weights = weights(l2_idx,:);
                    phase = phase(z_idx,:);
                    if isempty(noisevar)
                        obj.runLDscore;
                    end
                end
                
                % Normalizing z scores and noisevar
                obj.normalizer = obj.normalize;
                
            else
                error('input to DATA should be either a MODEL, the path to a directory, or a cell array of paths to files')
            end
            
        end
        
        function runLDscore(obj)
            LDSCout = CLDSC(obj);
            obj.sigmaEps = LDSCout.intercept;
            obj.ldsc=LDSCout;
        end
        
        function normalizer = normalize(obj)
            % covariance of alpha
                trait_cov=obj.z'*obj.z/length(obj.z)-obj.sigmaEps;
                
                % 1/sqrt(var(alpha))
                normalizer=1./sqrt(diag(trait_cov));
                
                % normalize sumstats so var(E(ZZ|alpha))=1
                obj.z=obj.z.*normalizer';
                
                % correlation of epsilon
                obj.sigmaEps=diag(normalizer)*obj.sigmaEps*diag(normalizer);
                
        end
        
        function obj = simulate(obj,model,mm)
            obj.beta = model.simulate(mm);
            obj.z = obj.beta + randn(mm,obj.noTraits)*chol(obj.sigmaEps);
            
        end
        
        function a = get.noTraits(obj)
            a = size(obj.z,2);
        end
        function a = get.noSNPs(obj)
            a = size(obj.z,1);
        end
        function a = get.noBlocks(obj)
            a = length(obj.blocks);
        end
        
        function a = get.cov(obj)
            a = cov(obj.z) - obj.sigmaEps;
        end
    end
end
