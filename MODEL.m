classdef MODEL < handle
    properties
        cpts
        traits {iscell} % trait names
        prediction % vector of marginal effect size predictions
        additive_cpts
        df
    end
    
    properties (Dependent)
        noCpts
        noIngs
        noTraits
        cov
        grad
    end
    
    methods
        function obj = MODEL(cpts,additive_cpts)
            if nargin==0
                return;
            end
            obj.cpts=cpts;
            if nargin < 2
                obj.additive_cpts = 1;
            else
                obj.additive_cpts = additive_cpts;
            end
            obj.df = 0;
        end
        
        %%  simulate
        function [beta, whichSNPs] = simulate(obj,mm)
            
            cptWeight = cellfun(@sum,{obj.cpts.ww});
            if all(cptWeight == 0)
                error('Weights must not all be zero to simulate from model')
            end
            
            beta=zeros(mm,obj.noTraits);
            
            if obj.additive_cpts
                for i = 1:obj.noCpts
                    inCpt = rand(mm,1) < sum(obj.cpts(i).ww);
                    beta(inCpt,:) = beta(inCpt,:) + ...
                        obj.cpts(i).simulate(sum(inCpt));
                    whichSNPs{i} = find(inCpt == 1, length(inCpt));
                end
            else
                
                if sum(obj.ww) ~= 1
                    error('Weights must sum to 1 to simulate from non-additive model')
                end
                whichCpt = randsample(1:obj.noCpts,mm,true,cptWeight/sum(cptWeight));
                for i = 1:obj.noCpts
                    beta(whichCpt==i,:)=obj.cpts(i).simulate(sum(whichCpt==i));
                    whichSNPs{i} = find(whichCpt == i, length(whichCpt))';
                end
            end
            
            
        end
        
        
        %% fit
        % perform gradient descent to fit parameters of model
        function [rss_traj] = fit(obj,ecf,varargin)
            obj.traits = ecf.traits;
            
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'ecf', @(obj)isa(obj,'ECF'));
            
            addParameter(p, 'convergence_param', 1e-6);
            addParameter(p, 'gdsteps', 100);
            addParameter(p, 'nonneg_weights', 1);
            addParameter(p, 'stepsize_param', 1/2);
            addParameter(p, 'print_stuff', true);
            parse(p,obj,ecf,varargin{:});
            convergence_param=p.Results.convergence_param;
            gdsteps=p.Results.gdsteps;
            nonneg_weights=p.Results.nonneg_weights;
            stepsize_param = p.Results.stepsize_param;
            if stepsize_param <= 0 || stepsize_param >= 1
                error('stepsize_param should be in (0,1)')
            end
            
            if numel(obj.theta) == 0
                obj.project(ecf,nonneg_weights);
                rss_traj(1) = obj.rss(ecf);
                return
            end
            currentPhi = obj.project(ecf,nonneg_weights);
            obj.cpts.gradCompute(ecf, currentPhi); % compute initial gradient
            converged = 0;
            stepsizeOld=100;
            
            % gradient descent
            warning('off','MATLAB:rankDeficientMatrix')
            if p.Results.print_stuff
                fprintf('Iteration: ')
            end
            counter=1;
            
            %             rss_traj = zeros(gdsteps,1);
            rss_traj(1) = obj.rss(ecf);
            
            while converged<2
                if p.Results.print_stuff
                    fprintf('%d ',counter)
                end
                counter=counter+1;
                if counter==gdsteps
                    if gdsteps > 10
                        warning('\n Failed to converge after %d steps',counter)
                    end
                    break;
                end
                
                rssOld = obj.rss(ecf);
                %gradOld=obj.grad;
                thetaOld = obj.theta;
                stepsize = stepsizeOld / stepsize_param;
                
                obj.theta(obj.theta - stepsize*obj.grad');
                
                % https://math.stackexchange.com/questions/373868/optimal-step-size-in-gradient-descent
                while obj.rss(ecf) > rssOld - stepsize_param * stepsize * sum(obj.grad.^2)
                    stepsize = stepsize * stepsize_param;
                    obj.theta(thetaOld - stepsize*obj.grad');
                end
                stepsizeOld = stepsize;
                
                % only call project again after finding good stepsize
                xbeta = obj.project(ecf,nonneg_weights);
                %disp(rssOld - obj.rss(ecf))
                % compute new gradient
                dX=obj.cpts.gradCompute(ecf, xbeta);
                
                
                if isnan(stepsize)
                    fprintf('\n')
                    warning('Gradient descent became stuck on iteration %d',counter)
                    break;
                end
                
                rssCur = obj.rss(ecf);
                converged = converged + (abs(rssOld-rssCur) / rssCur < convergence_param);
                
                rss_traj(counter) = obj.rss(ecf);
            end
            if p.Results.print_stuff
                fprintf('%d\n', counter)
            end
            %             rss_traj = rss_traj(1:counter);
            
            % df after calling fit
            obj.df=min(obj.df + numel(obj.theta),...
                numel(ecf.phiP));
            
        end
        
        %% project
        % modelPhi is either a vector (mixture case) or a matrix (additive
        % case) corresponding to the characteristic function of the model
        % (mixture case) evaluated at ecf.T, or the characteristic function
        % of each component (additive case).
        % no_iter_project is number of iterations for iterative projection
        % caluclation (OK if it doesn't go all the way to convergence
        % within each step of gradient descent)
        function xbeta = project(obj,ecf,nonneg_weights,sum_weights_one,no_iter_project)
            warning('off','MATLAB:singularMatrix')
            warning('off','MATLAB:illConditionedMatrix')
            warning('off','MATLAB:rankDeficientMatrix')
            if nargin<3
                nonneg_weights = 1;
            end
            if nargin < 4
                sum_weights_one = 1;
            end
            if nargin < 5
                no_iter_project = 5;
            end
            
            % additive model
            if obj.additive_cpts
                xc = cell(obj.noCpts,1);
                beta = xc;
                if ~isempty([obj.cpts.A])
                    warning('Constraints are ignored in additive case')
                end
                
                % phi_epsilon(t)
                yhat = exp(-1/2 * sum(ecf.T.*(ecf.sigmaEps * ecf.T),1)' ) ;
                
                for cc = 1:obj.noCpts
                    % phi_k(t) - 1 (no noise factor)
                    xc{cc} = obj.cpts(cc).phi(ecf, zeros(obj.noTraits));
                    
                    beta{cc} = obj.cpts(cc).ww;
                    
                    % yhat = phi_epsilon * ((phi_1 - 1)*ww_1 + 1) * ...
                    yhat = yhat .* ( xc{cc} * beta{cc} + 1);
                end
                
                
                for rep = 1:no_iter_project
                    for cc = 1:obj.noCpts
                        if sum_weights_one
                            AA = ones(1,obj.cpts(cc).noIngs);
                            conval = 1;
                        else
                            AA = [];
                            conval = [];
                        end
                        %current_err = sum((ecf.phiP + 1 - ecf.P * yhat).^2);
                        yhat = yhat ./ (xc{cc} * beta{cc} + 1 );
                        beta{cc} = constrained_regression(...
                            ecf.P * (xc{cc} .* yhat) , ...
                            ecf.phiP - ecf.P * (yhat - 1),...
                            ones(size(ecf.phiP)),nonneg_weights,[],AA,conval);
                        yhat = yhat .* (xc{cc} * beta{cc} + 1);
                    end
                end
                
                xbeta = cellfun(@(x,b){1 + x*b},xc,beta);
                xbeta = [xbeta{:}];
                
                beta = vertcat(beta{:});
                % mixture model
            else
                xr=ecf.P * obj.cpts.phi(ecf);
                beta = constrained_regression(xr, ecf.phiP,...
                    ones(size(ecf.phiP)),nonneg_weights,[],[obj.cpts.A],[]);
                xbeta = 1 + obj.cpts.phi(ecf)*beta;
            end
            
            obj.cpts.setww(beta);
            
            % no. dimensions
            if numel(obj.ww) - rank([obj.cpts.A]) > numel(ecf.phiP)
                warning('Model has too many parameters compared with the number of sampling-time combinations')
            end
            
            if isempty(obj.traits)
                obj.traits = ecf.traits;
            end
            % df after calling project: # observations - whichever smaller:
            % (a) number of weights minus number of constraints on those
            % weights, ie the maximum column rank
            % (b) the number of sampling times, ie the row rank
            obj.df= min(numel(obj.ww) - rank([obj.cpts.A]), numel(ecf.phiP));
            
        end
        
        %% initialize_fit
        % Function to iterate through different initializations of theta
        % and call fit - has different ways of initializing theta,
        % including starting at the initial value of theta (can set it to
        % be the true theta for simulations), randomized theta, a mixture
        % of true and random initializations, and initializations near the
        % true theta
        function [rss_traj,inittime,gdtime] = initialize_fit(obj, ecf, varargin)
            
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'ecf', @(obj)isa(obj, 'ECF'));
            
            addParameter(p, 'niter', 100); % number of initializations of gradient descent
            addParameter(p, 'nrot', 5); % Max number of rotations
            addParameter(p, 'gdsteps', 100); % Max number of gradient descent steps
            addParameter(p, 'init_gdsteps', 3); % Max number of gradient descent steps during initialization
            addParameter(p, 'nonneg', 1');% Constrain model to have nonnegative weights
            addParameter(p, 'method', 'covmat'); % Initialization method
            addParameter(p, 'max_offset', 1); % goes with method "near_orig"; how close to initial param values to initialize
            addParameter(p, 'fix_cpt', []);% indices of components for which to fix the params at their initial value
            addParameter(p, 'dampening_factor', 0);% Dampens changes in estimated residuals-covariance matrix to improve convergence
            addParameter(p, 'orig_subset', []);% goes with method "orig_subset"; indices of cpts to fix params at initial value; all other cpts are randomized
            addParameter(p, 'datacov', []);% covariance matrix used with 'covmat' initialization method
            addParameter(p, 'nsamples', []);% goes with datacov option + 'covmat'
            addParameter(p, 'convergence_param', 1e-6);% threshold for convergence
            addParameter(p, 'rotation_tolerance', 1e-2); % ECF rotation parameter
            addParameter(p, 'print_stuff', true);
            
            parse(p,obj,ecf,varargin{:});
            
            niter=p.Results.niter;
            nrot=p.Results.nrot;
            gdsteps=p.Results.gdsteps;
            init_gdsteps=p.Results.init_gdsteps;
            nonneg=p.Results.nonneg;
            method=p.Results.method;
            max_offset=p.Results.max_offset;
            fix_cpt=p.Results.fix_cpt;
            orig_subset=p.Results.orig_subset;
            datacov=p.Results.datacov;
            nsamples=p.Results.nsamples;
            convergence_param=p.Results.convergence_param;
            
            methodlist = ["orig", "orig_rand", "random", "near_orig", "orig_subset", "covmat"];
            if ~any(strcmp(methodlist,method))
                error("Error: Need to method must be one of the following: 'orig', 'random', 'orig_rand', 'near_orig', 'orig_subset', 'covmat'");
            end
            
            % copying original components
            components = copy(obj.cpts);
            truetheta = obj.theta;
            
            % iterate until finding a good fit
            best=inf;
            
            tic;
            for i=1:niter
                % To do: rewrite long if-else
                switch method
                    case "orig"
                        % Initialize theta to its initial value
                        obj.cpts = components;
                    case "orig_rand"
                        % On first iteration, initialize theta to initial value
                        if i == 1
                            obj.cpts = components;
                        else
                            if isempty(datacov)
                                error("Need to supply data covariance for covmat method");
                            end
                            if isempty(nsamples)
                                nsamples = size(datacov,1);
                            end
                            for cc=1:obj.noCpts
                                cpt = obj.cpts(cc);
                                if cpt.name == "rank_one" || cpt.name == "full_rank"
                                    cpt.settheta(sample_datacov(datacov, cpt.name,'nsamples',nsamples));
                                else
                                    cpt.randomize;
                                end
                            end
                        end
                    case "random"
                        % Randomize theta on every iteration
                        obj.cpts.randomize;
                    case "near_orig"
                        % Initialize theta near the initial value
                        newtheta = truetheta + rand(size(truetheta))*2*max_offset - max_offset;
                        obj.theta(newtheta);
                    case "orig_subset"
                        if isempty(orig_subset)
                            warning('Did not pass in any components to use initial vals, so randomizing all components.');
                        end
                        obj.cpts(orig_subset) = copy(components(orig_subset));
                        obj.cpts(setdiff(1:obj.noCpts,orig_subset)).randomize;
                    case "covmat"
                        % initializes pleiotropic components to be drawn from
                        % multivariate normal with covariance from data
                        if isempty(datacov)
                            error("Need to supply data covariance for covmat method");
                        end
                        if isempty(nsamples)
                            nsamples = size(datacov,1);
                        end
                        for cc=1:obj.noCpts
                            cpt = obj.cpts(cc);
                            if cpt.name == "rank_one" || cpt.name == "full_rank"
                                cpt.settheta(sample_datacov(datacov, cpt.name,'nsamples',nsamples));
                            else
                                cpt.randomize;
                            end
                        end
                end
                
                
                if ~isempty(fix_cpt)
                    obj.cpts(fix_cpt) = copy(components(fix_cpt));
                end
                
                obj.fit(ecf,'nonneg_weights',nonneg,'gdsteps',init_gdsteps,...
                    'convergence_param',convergence_param,'print_stuff',p.Results.print_stuff);
                rssCur = obj.rss(ecf);
                if rssCur < best
                    besttheta=obj.theta;
                    best=rssCur;
                    bestcpts = copy(obj.cpts);
                end
                
                if method == "orig"
                    break
                end
                
            end
            inittime=toc;
            
            obj.cpts = bestcpts;
            
            % Rotating and fitting until the solution converges or until
            % max rotations is reached
            oldcov=obj.cov;
            covdiffthresh=0.1;
            
            tic;
            for ii=1:nrot
                obj.rotateECF(ecf,'dampening_factor',p.Results.dampening_factor,...
                    'tol', p.Results.rotation_tolerance);
                rss_traj = obj.fit(ecf,'gdsteps',gdsteps,...
                    'convergence_param',convergence_param,...
                    'print_stuff',p.Results.print_stuff);
                
                
                if ~isempty(fix_cpt)
                    obj.cpts(fix_cpt) = copy(components(fix_cpt));
                end
                
                cdiff = covdiff(oldcov,obj.cov);
                oldcov = obj.cov;
                if cdiff <= covdiffthresh
                    break
                end
            end
            gdtime=toc;
            
            if ii==nrot
                warning('Estimates did not converge');
            end
            
        end
        
        %% rotateECF
        % Computes new rotation/projection (ecf.P) based on the model
        function [phi, phicombos, cpt_phi] = rotateECF(obj, ecf, varargin)
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'ecf', @(obj)isa(obj, 'ECF'));
            
            addParameter(p, 'tol', 0);
            addParameter(p, 'fudge_factor', 0);
            addParameter(p, 'dampening_factor', 0.5);
            
            parse(p,obj,ecf,varargin{:});
            
            tol=p.Results.tol;
            fudge_factor=p.Results.fudge_factor;
            dampening_factor=p.Results.dampening_factor;
            % pairwise sums + differences of sampling times
            Tsums = zeros(obj.noTraits,ecf.noT*(ecf.noT+1)/2);
            Tdiffs = zeros(obj.noTraits,ecf.noT*(ecf.noT+1)/2);
            for kk=1:obj.noTraits
                Tsums(kk,:) = triuind(ecf.T(kk,:) + ecf.T(kk,:)');
                Tdiffs(kk,:) = triuind(ecf.T(kk,:) - ecf.T(kk,:)');
            end
            
            % phi(t1+t2) + phi(t1-t2)
            
            sigmaEps = ecf.sigmaEps + fudge_factor*eye(length(ecf.sigmaEps));
            
            if obj.additive_cpts
                phicombos =  ones(size(Tsums,2),1);
                phi = ones(ecf.noT,1);
                for cc = 1 : obj.noCpts
                    cpt_phi = (obj.cpts(cc).phi(Tsums, sigmaEps) ...
                        + obj.cpts(cc).phi(Tdiffs, sigmaEps))/2;
                    phicombos = phicombos .* (cpt_phi*obj.cpts(cc).ww + 1);
                    phi = phi .* (obj.cpts(cc).phi(ecf)*obj.cpts(cc).ww + 1);
                end
            else
                cpt_phi = (obj.cpts.phi(Tsums, sigmaEps) + ...
                    obj.cpts.phi(Tdiffs, sigmaEps))/2;
                phicombos = cpt_phi*obj.ww + 1;
                phi = obj.cpts.phi(ecf) * obj.ww + 1;
            end
            
            % 1/2 * (phi(t1+t2) + phi(t1-t2)) - phi(t1)phi(t2)
            
            phicombos = phicombos - triuind(phi * phi');
            
            % convert back to matrix
            phicombos = triuind(phicombos);
            phicombos = phicombos + phicombos' - diag(diag(phicombos));
            
            
            % calculate rotation
            if dampening_factor > 0
                oldcombos = ecf.P' * diag(1./sum(ecf.P'.^2).^2) * ecf.P;
                phicombos = (1 - dampening_factor) * phicombos + ...
                    dampening_factor * oldcombos;
                phicombos = (phicombos + phicombos') / 2;
            end
            
            [U,S] = eig(phicombos);
            incl = diag(S) > tol;
            Sinv = 1./sqrt(diag(S));
            ecf.P = diag(Sinv(incl)) * U(:,incl)';
            
        end
        
        %% test_cpts
        % Function to test the significance of each component by removing
        % it, or replacing uncorrelated components with trait_specific;
        % also tests the significance of the nullmodels against any other
        % model
        function [pval, pval_general_H1] = test_cpts(obj, ecf, varargin)
            
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'ecf', @(obj)isa(obj, 'ECF'));
            addParameter(p, 'traitIdx', 1:obj.noTraits);
            addParameter(p, 'nrot', 1);
            addParameter(p, 'gdsteps', 1e3);
            addParameter(p, 'whichCpts', 1:obj.noTraits);
            addParameter(p, 'tol', 0.005);
            parse(p,obj,ecf,varargin{:});
            
            nrot=p.Results.nrot;
            gdsteps=p.Results.gdsteps;
            whichCpts=p.Results.whichCpts;
            
            % First testing fitted model against the general model
            pval_general_H1(obj.noCpts+1) = obj.test(ecf);
            
            % Whether to replace or remove uncorrelated components
            orig_cpts = copy(obj.cpts);
            for kk=whichCpts
                nullcpts = copy(orig_cpts([1:kk-1, kk+1:end]));
                
                nullmodel = MODEL(nullcpts,obj.additive_cpts);
                
                nullmodel.fit(ecf,'gdsteps',gdsteps);
                for r = 1:nrot
                    nullmodel.rotateECF(ecf,'tol',p.Results.tol);
                    nullmodel.fit(ecf,'gdsteps',gdsteps);
                end
                
                
                % Fitting full model based on rotation on nullmodel
                obj.fit(ecf,'gdsteps',gdsteps);
                
                pval(kk)=nullmodel.test(obj,'ecf',ecf);
                pval_general_H1(kk) = nullmodel.test(ecf);
                
                obj.cpts = copy(orig_cpts);
                
            end
            
        end
        
        %%
        % test model as null vs some alternative: either another model, or
        % an ecf object
        function [pval, Fstat, convergence_tol] = test(H0,input,varargin)
            
            p=inputParser;
            
            addRequired(p, 'H0', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'input', @(obj)isa(obj,'MODEL') || isa(obj,'ECF'));
            
            addParameter(p, 'ecf', []);
            
            parse(p,H0,input,varargin{:});
            
            ecf = p.Results.ecf;
            
            if isa(input, 'MODEL')
                if isempty(ecf)
                    error("If inputting a model, must also input ecf with 'ecf' flag");
                end
                noData = numel(ecf.phiP);
            else
                ecf = input;
                noData = numel(ecf.phiPBlocks);
            end
            
            % test vs alternative model
            if isa(input,'MODEL')
                H1=input;
                if H0.df >= H1.df
                    warning('Attempting to compare a null that is the same size or larger than the alternative will produce NaN p-values')
                    Fstat = 0;
                    pval = NaN;
                    return;
                end
                
                H0_rss = H0.rss(ecf);
                H1_rss = H1.rss(ecf);
                
                if H0_rss < H1_rss
                    warning('Null model has lower RSS than alternative');
                end
                
                Fstat=(H0_rss-H1_rss)/(H1.df-H0.df) / ...
                    (H1_rss / (noData - H1.df));%%%
                pval=fcdf(Fstat,H1.df-H0.df,noData - H1.df,'upper');
                
                H1_df = H1.df;
                
                % test vs generalized alternative
            elseif isa(input,'ECF')
                
                H1_rss = sum(sum((input.phiPBlocks - input.phiP).^2));
                H1_df = numel(input.phiP);
                if  H0.df >= H1_df
                    warning('Attempting to compare a null that is the same size or larger than the alternative will produce NaN p-values')
                    Fstat = 0;
                    pval = NaN;
                    return;
                end
                
                H0_rss = H0.rss(ecf,1);
                
                Fstat=(H0_rss-H1_rss)/(H1_df - H0.df) / ...
                    (H1_rss / (noData - H1_df));%%%
                pval=fcdf(Fstat,H1_df - H0.df, noData - H1_df,'upper');
                
            else
                error('Input to MODEL.test should be either a MODEL or an ECF')
            end
            
            convergence_tol=0.1*finv(.001,H1_df-H0.df,noData-H1_df)*(H1_df-H0.df)/(noData-H1_df);
        end
        
        %%
        % Calculate expected value of marginal effect sizes
        % Outputs:
        % - alpha: predicted effect sizes
        % - MAP_scalars: scalars of the components based on MAP ingredient assignment
        % - posterior_scalars: posterior variance calculated from the
        %     posterior likelihood
        function [alpha, MAP_scalars, posterior_scalars] = predict(obj,data,varargin)
            
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            addRequired(p, 'data', @(obj)isa(obj, 'DATA'));
            
            addParameter(p, 'traitIdx', 1:data.noTraits); % traits to predict effect sizes for
            addParameter(p, 'whichSNPs', 1:data.noSNPs); % SNPs to predict effect sizes for
            addParameter(p, 'marg', 0); % if 1, generate predictions based on each trait's univariate distribution
            addParameter(p, 'minWeight', 0); % weight product threshold to select ingredient combinations
            addParameter(p, 'sigma_extra', zeros(data.noTraits)); % amount of noise to transfer from noisevar to the covariance matrix
            
            parse(p,obj,data,varargin{:});
            
            traitIdx=p.Results.traitIdx;
            marg=p.Results.marg;
            whichSNPs = p.Results.whichSNPs;
            minWeight = p.Results.minWeight;
            sigma_extra = p.Results.sigma_extra;
            
            x = data.z(whichSNPs,traitIdx);
            sigmaEps = data.sigmaEps(traitIdx,traitIdx);
            
            weightproducts = obj.cpts(1).ww;
            for ii = 2:obj.noCpts
                weightproducts = weightproducts * obj.cpts(ii).ww';
                weightproducts = weightproducts(:);
            end
            
            noChoices = [obj.cpts.noIngs];
            linInds = find(weightproducts > minWeight);
            weightproducts = weightproducts(weightproducts > minWeight);
            indArray = indexFn(linInds-1,noChoices)+1;
            
            if ~marg
                xpm = zeros(size(x));
                likelihood = zeros(size(xpm,1),1);
            else
                xpm = zeros(size(x));
                likelihood = xpm;
            end
            
            cpt_sigma = {obj.cpts.S};
            cpt_scalars = {obj.cpts.scalars};
            
            %             posterior_likelihood = zeros(size(x,1),length(indArray));
            max_pl = ones(size(x,1),1)*-Inf;
            indMax = ones(size(x,1),1)*-Inf;
            for ii = 1:length(indArray)
                kk = num2cell(indArray(ii,:));
                sigma_cells = cellfun(@(sigma,scalar,kk)sigma * scalar(kk),...
                    cpt_sigma, cpt_scalars, kk,'UniformOutput',0);
                sigma_sum = sum(cat(3,sigma_cells{:}),3);
                
                if ~marg
                    [xpm_temp, likelihood_temp] = normpm(x,sigma_sum + sigma_extra,...
                        data.sigmaEps - sigma_extra);
                    
                    xpm = xpm + likelihood_temp .* weightproducts(ii) .* xpm_temp';
                    likelihood = likelihood + likelihood_temp * weightproducts(ii);
                    
                    posterior_likelihood{ii} = likelihood_temp * weightproducts(ii);
                    [max_pl_temp, ind_max_pl_temp] = max([max_pl posterior_likelihood{ii}],[],2);
                    max_pl = max_pl_temp;
                    indMax(ind_max_pl_temp == 2) = ii;
                else
                    for t=1:obj.noTraits
                        [xpm_temp, likelihood_temp] = normpm(x(:,t),sigma_sum(t,t) + sigma_extra(t,t),...
                            data.sigmaEps(t,t) - sigma_extra(t,t));
                        
                        xpm(:,t) = xpm(:,t) + likelihood_temp .* weightproducts(ii) .* xpm_temp';
                        likelihood(:,t) = likelihood(:,t) + likelihood_temp * weightproducts(ii);
                        
                        posterior_likelihood{ii} = likelihood_temp * weightproducts(ii);
                        [max_pl_temp, ind_max_pl_temp] = max([max_pl posterior_likelihood{ii}],[], 2);
                        max_pl = max_pl_temp;
                        indMax(ind_max_pl_temp == 2) = ii;
                        
                    end
                    
                end
                
                
            end
            %             [~, indMax ] = max(posterior_likelihood');
            MAP_assignment = indArray(indMax,:);
            for cpt = 1:obj.noCpts
                MAP_scalars(:,cpt) = obj.cpts(cpt).scalars(MAP_assignment(:,cpt));
            end
            posterior_likelihood_sum = zeros(size(x,1), obj.noCpts, max([obj.cpts.noIngs]));
            posterior_scalars = zeros(size(x,1),obj.noCpts);
            for ii = 1:obj.noCpts
                for jj=1:obj.cpts(ii).noIngs
                    if sum(indArray(:,ii) == jj) > 0
                        temp_idcs = find(indArray(:,ii) == jj);
                        for kk = 1:length(temp_idcs)
                            posterior_likelihood_sum(:,ii,jj) = posterior_likelihood_sum(:,ii,jj) + posterior_likelihood{temp_idcs(kk)};
                        end
                    end
                    
                end
                temp = reshape(posterior_likelihood_sum(:,ii,:),size(x,1),max([obj.cpts.noIngs]));
                posterior_scalars(:,ii) = ...
                    temp./sum(temp,2) * ...
                    obj.cpts(ii).scalars';
            end
            
            obj.prediction = xpm ./ likelihood;
            alpha=obj.prediction;
            
        end
        
        % Function to compute SSE; if a data object is input, it also
        % computes the predictions, otherwise a second set of effect sizes
        % needs to be provided
        function sse = sse(obj, alpha1, alpha2)
            if nargin < 3
                alpha2 = obj.prediction;
            end
            
            sse = sum((alpha1-alpha2).^2);
        end
        
        % set obj.cpts.theta to specified values; also returns it
        function th = theta(obj,newTheta)
            if nargin > 1
                obj.cpts.settheta(newTheta);
            end
            th = vertcat(obj.cpts.theta);
            
        end
        
        % set obj.cpts.ww to specified values; also returns it
        function w = ww(obj,newWw)
            if nargin > 1
                obj.cpts.setww(newWw);
            end
            w = vertcat(obj.cpts.ww);
            
        end
        
        % function to get rss
        function rss = rss(obj, ecf, useBlocks)
            if nargin < 3
                useBlocks = 0;
            end
            % phi_epsilon(t)
            modelPhi = exp(-1/2*sum(ecf.T.*(ecf.sigmaEps*ecf.T),1)');
            
            if obj.additive_cpts
                for cc = 1:obj.noCpts
                    xx = obj.cpts(cc).phi(ecf, 0);
                    beta = obj.cpts(cc).ww;
                    modelPhi = modelPhi .* (1 + xx*beta);
                end
            else
                xx = obj.cpts.phi(ecf, 0);
                beta = obj.ww;
                modelPhi = modelPhi .* (1 + xx*beta);
            end
            
            if useBlocks
                rss=sum(sum( (ecf.phiPBlocks - ecf.P * (modelPhi-1) ).^2));
            else
                rss=sum(sum( (ecf.phiP - ecf.P * (modelPhi-1) ).^2));
            end
        end
        
        function h = heatmap(obj,cpt,tnames,titles,plot,cmap)
            if nargin < 5
                plot = 1;
                cmap = 1;
            elseif nargin < 6
                cmap = 1;
            end
            
            D = diag(1./sqrt(diag(obj.cov)));
            h = D*obj.cpts(cpt).cov*D;
            h = round(h,2);
            
            if plot
                h1 = heatmap(tnames,tnames,h);
                h1.FontSize = 14;
                colormap(bluewhitered(256))
                if ~cmap
                    colorbar off
                end
                caxis([-1 1])
                title(titles)
            end
        end
        
        function [ax] = h2plot(obj,varargin)
            p=inputParser;
            
            addRequired(p, 'obj', @(obj)isa(obj,'MODEL'));
            
            addParameter(p, 'traitNames', obj.traits,@iscell);
            addParameter(p, 'whichCpts', 1:obj.noCpts, @(x)size(x,1)==1);
            addParameter(p, 'cptNames', []);
            addParameter(p, 'whichTraits', 1:obj.noTraits);
            
            parse(p,obj,varargin{:});
            
            if isempty(p.Results.cptNames)
                %                 cptNames = {obj.cpts(p.Results.whichCpts).name};
                cptnums = string(1:obj.noCpts);
                cptNames = cellstr(arrayfun(@(x)join(['cpt',x]),cptnums));
            else
                cptNames = p.Results.cptNames;
            end
            
            for i=1:obj.noCpts
                cpth2(:,i)=diag(obj.cpts(i).cov);
            end
            cpth2=cpth2(p.Results.whichTraits,:);
            cpth2 = cpth2./sum(cpth2,2);
            cpth2 = round(cpth2,2);
            heatmap(cptNames,p.Results.traitNames(p.Results.whichTraits),cpth2(:,p.Results.whichCpts))
            ax=gca;
            colormap(bluewhitered(256));
            caxis([-1 1]);
        end
        
        function [tables,fields] = print_cpt_data(obj)
            fields = ["theta","ww","cov","s"];
            for ii = 1:length(fields)
                varNames = arrayfun(@(x)join(['cpt',num2str(x)]),1:obj.noCpts,'UniformOutput',false);
                switch ii
                    case 1
                        datasz = max(cellfun(@(x)length(x),{obj.cpts.theta}));
                    case 2
                        datasz = max(cellfun(@(x)length(x),{obj.cpts.ww}));
                    case 3
                        datasz = max(cellfun(@(x)length(triuind(x)),{obj.cpts.cov}));
                    case 4
                        datasz = max(cellfun(@(x)length(triuind(x)),{obj.cpts.S}));
                end
                sz = [datasz obj.noCpts];
                varTypes = repmat("double",1,obj.noCpts);
                temp = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
                for cc = 1:obj.noCpts
                    switch ii
                        case 1
                            vals = obj.cpts(cc).theta;
                            if length(vals) > 1
                                temp(1:length(vals),cc) = num2cell(vals);
                            end
                        case 2
                            vals = obj.cpts(cc).ww;
                            temp(1:length(vals),cc) = num2cell(vals);
                        case 3
                            vals = triuind(obj.cpts(cc).cov);
                            temp(1:length(vals),cc) = num2cell(vals);
                        case 4
                            vals = triuind(obj.cpts(cc).S);
                            temp(1:length(vals),cc) = num2cell(vals);
                    end
                end
                tables{ii} = temp;
            end
        end
        
        % total covariance of model
        function a = get.cov(obj)
            a = sum(cat(3,obj.cpts.cov),3);
        end
        
        function a = get.noCpts(obj)
            a = numel(obj.cpts);
        end
        function a = get.noIngs(obj)
            a = sum([obj.cpts.noIngs]);
        end
        % calculate gradient of each component and project onto tangent
        % space of constraint level sets
        function a = get.grad(obj)
            a = horzcat(obj.cpts.grad);
            dCon = horzcat(obj.cpts.dA);
            % residualize gradient on the gradient of the constraints
            if ~isempty(dCon)
                a = a - (a/dCon)*dCon;
            end
        end
        
        function a = get.noTraits(obj)
            if obj.noCpts > 0
                nT=[obj.cpts.noTraits];
                a = nT(1);
                if any(nT ~= nT(1))
                    error('All components should have the same number of traits')
                end
            end
        end
    end
end
