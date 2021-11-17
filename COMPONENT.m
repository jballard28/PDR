classdef COMPONENT < matlab.mixin.Copyable
    properties
        theta % column vector of parameters for this component
        scalars % one scalar for each ingredient
        name % name of component (just for convenience)
        ww % weight of each ingredient
        thetaCon % upper and lower bounds on theta
        noTraits
        grad
        conType
        fixedParam
    end
    
    properties (Dependent)
        A % matrix of constraint values, no. constraints x no. ingredients
        dA % gradient of A
        noIngs % length(scalars)
        noCons % number of coefficient constraints
        noParams
        cov % S*dot(scalars,ww)
        S % covariance matrix of the component
    end
    
    properties (Access = private)
        cptfn % used to compute gradients
        vCon % used to compute constraint coefficients
        type
    end
    
    properties (Dependent, Access = private)
        F
        dF
    end
    
    methods
        function obj = COMPONENT(type,no_traits,scalars,fixedParam,varargin)
            if nargin==0
                return;
            end
            p=inputParser;
            
            % which type of component (member of cpts folder)
            addRequired(p, 'type', @isstr);
            
            % scalars, one per ingredient
            addRequired(p, 'no_traits', @(x)isnumeric(x) & isscalar(x));
            
            % scalars, one per ingredient
            addRequired(p, 'scalars', @(x)isnumeric(x) & isvector(x) & size(x,1)==1);
            
            % parameter (different from theta) that isn't learned
            addRequired(p, 'fixedParam', @(x)true)
            
            % parameters; if unspecified, these are randomly initialized
            [~, ~, randomParams, obj.thetaCon] = component_function(type,no_traits,fixedParam,1);
            
            addOptional(p, 'theta', randomParams, ...
                @(x)isnumeric(x) & all(size(x)==size(randomParams)));
            
            % weight of each ingredient
            addOptional(p, 'ww', zeros(size(scalars))', ...
                @(x)isnumeric(x) & size(x,1) == size(scalars,2));
            
            
            % name of component (just for convenience)
            addOptional(p, 'name', type, @isstr)
            
            % option to generate heritability constraints. Options allowed:
            % 'none' (no constraint), 'diagonal' (constrain h2, but not
            % genetic correlation), 'full' (constrain h2 + rg)
            addOptional(p, 'h2ConOpt', 'none', @isstr)
            
            
            % custom constraints that can be used to constrain how much
            % heritability/genetic covariance is explained by different
            % components. First two columns (interchangeable) are indices
            % of a trait pair (entry of gencov matrix), with duplicates
            % allowed; third column is a coefficient.
            % Not yet implemented.
            %addOptional(p, 'customCon', zeros(0,3), @(x)size(x,2)==3)
            
            parse(p,type,no_traits,scalars,randomParams,varargin{:});
            
            obj.theta=p.Results.theta;
            obj.noTraits=no_traits;
            obj.type = type;
            obj.conType = p.Results.h2ConOpt;
            obj.fixedParam = fixedParam;
            
            % f(t) is equal to t'*Sigma*t, where Sigma is a function of
            % theta. df(t) is the gradient of this wrt theta.
            obj.cptfn = @(th)component_function(type, no_traits, fixedParam, th);
            
            % constraint function of component
            i=strcmp(p.Results.h2ConOpt,{'none','diagonal','full'});
            
            if i(1)
                obj.vCon=zeros(no_traits,0);
            elseif i(2)
                obj.vCon=eye(no_traits);
            elseif i(3)
                error('Not yet implemented')
            else
                error('h2ConOpt must be either none, diagonal, or full')
            end
            
            %obj.Afn = @(S)[ones(no_traits,1).*p.Results.weightSum, S.*p.Results.h2sum];
            
            obj.scalars = scalars;
            obj.ww = p.Results.ww;
            obj.name = p.Results.name;
            
        end
        
        function newCon(obj,v)
            obj.vCon = [obj.vCon;v];
        end
        
        function delCon(obj,kk)
            if nargin < 1
                kk = 1:obj.noCons;
            end
            obj.vCon = obj.vCon(setdiff(1:end,kk),:);
        end
        
        % set theta to a new value while respecting parameter constraints
        function settheta(obj,thetaNew)
            counter=0;
            if numel(thetaNew)~=length(vertcat(obj.theta))
                error('new theta has wrong length for assignment')
            end
            for i = 1:numel(obj)
                theta_i=thetaNew(counter+1:counter+numel(obj(i).theta));
                counter=counter+numel(obj(i).theta);
                if numel(obj(i).thetaCon) > 0
                    obj(i).theta=min(obj(i).thetaCon(:,2),max(obj(i).thetaCon(:,1),theta_i));
                else
                    obj(i).theta = theta_i;
                end
            end
        end
        
        % evaluates the gradient of the RSS wrt theta
        function dX = gradCompute(obj,ecf,xbeta,additive_cpts)
            if nargin < 4
                additive_cpts = size(xbeta,2) > 1;
            end
            if additive_cpts
                if size(xbeta,2) ~= numel(obj)
                    error('Number of components should equal number of columns of xbeta')
                end
            end
            if size(xbeta,1) ~= ecf.noT
                error('xbeta should have a row for each sampling time')
            end
            
            % diagonal elements of T'*sigmaEps*T
            noisefactor=exp(-1/2*sum(ecf.T.*(ecf.sigmaEps*ecf.T),1)'); 
            
            % current value of the characteristic fn evaluated at ecf.T
            xbetaProd = prod(xbeta,2) .* noisefactor;
            
            for i=1:numel(obj)
                % evaluate the gradient of T'*S(theta)*T
                dFofT=obj(i).dF(ecf.T); % size(T,2) x noParams
                FofT=obj(i).F(ecf.T); % size(T,2) x 1
                
                % iterate over ingredients to compute gradient of
                % phi_theta(T)*ww (OK 4-8-21)
                dX=zeros(size(ecf.T,2),obj(i).noParams);
                for kk=1:obj(i).noIngs
                    dX = dX + obj(i).ww(kk) * ...
                        exp( -1/2 * (obj(i).scalars(kk)*FofT ) )... %deriv outside eval @ inside
                        .* ( -1/2 * obj(i).scalars(kk) * dFofT );%deriv inside
                end
                
                % additive case: dX gets multiplied by the factors from the
                % other components
                if additive_cpts
                    dX = dX .* xbetaProd ./ xbeta(:,i);
                else
                    dX = dX .* noisefactor;
                end
                
                % 
                
                % multiply gradient by projection matrix
                dX = ecf.P * dX;
                
                % gradient of the objective function
                obj(i).grad = -2 * (ecf.phiP - ecf.P * (xbetaProd - 1))'*dX;
            end
        end
        
        
        function randomize(obj)
            for i=1:numel(obj)
                [~, ~, randomParams] = component_function(obj(i).type,obj(i).noTraits,obj(i).fixedParam);
                obj(i).settheta(randomParams);
            end
        end
        
        function b=simulate(obj,mm)
            if numel(obj) ~= 1
                error('simulate method does not support nonscalar objects')
            end
            if mm > 0
                if any(obj.ww < 0)
                    error('Mixture weights should be nonnegative to simulate')
                end
                
                if all(obj.ww == 0)
                    error('Mixture weights should sum to a positive value')
                end
                
                % sample from MVN with covariance obj.S
                c = mvnrnd(zeros(1,obj.noTraits),obj.S,mm);
                
                % multiply by scalars chosen at random with weights proportional to ww
                if numel(obj.scalars) > 1
                    sc = randsample(sqrt(obj.scalars)',mm,true,obj.ww/sum(obj.ww));
                else
                    sc = sqrt(obj.scalars);
                end
                b=c.*sc;
                
            else
                b=zeros(0,obj.noTraits);
            end
            
        end
        
        % take a step down the gradient (supports nonscalar objects)
        function step(obj, stepsize)
            for i = 1:numel(obj)
                obj(i).settheta(obj(i).theta - stepsize * obj(i).grad');
            end
        end
        
        % set ww to a new value (supports nonscalar objects)
        function setww(obj, wwNew)
            if numel(wwNew) ~= sum([obj.noIngs])
                error('New weights vector has the wrong size')
            end
            counter = 0;
            for i=1:numel(obj)
                n = obj(i).noIngs;
                obj(i).ww(:) = wwNew(counter+1:counter+n);
                counter=counter+n;
            end
        end
        
        function a = phi(obj,input,noisevar)
            if isa(input,'ECF')
                samplingTimes = input.T;
                if nargin < 3
                    noisevar = input.sigmaEps;
                end
            elseif isnumeric(input)
                samplingTimes = input;
                if nargin < 3
                    noisevar = 0;
                end
            end
            
            % uses cellfun to support component arrays
            ingredient_quadratics=cellfun(@(a,f){a.*f(samplingTimes)},{obj.scalars},{obj.F});
            
            a = (exp( -1/2*( horzcat(ingredient_quadratics{:})... % like [scalars*t'*S(theta)*t]
                + sum( samplingTimes.*(noisevar * samplingTimes), 1 )' ) ) - 1); % like t'*sEps*t
        end
        
        % Calculate the mean of the posterior for a given ingredient
        % E(alpha|alpha_hat,G=g)
        function mean_posterior = mean_posterior(obj,x,sigmaEps,ing,traitIdx,marg)
            
            if ~marg
                det_thresh = 1e-14;
                s = obj.S*obj.scalars(ing);
                if (abs(det(s)) <= det_thresh)
                    scale = 1e-12;
                    s = s + eye(obj.noTraits)*scale;
                end
                
                temp = inv(s)+inv(sigmaEps);
                f=@(x)(temp\(sigmaEps\x))';
                a=cellfun(f,num2cell(x',1),'UniformOutput',false);
                mean_posterior = vertcat(a{:});
                
            else
                % For a given ingredient
                mean_posterior = zeros(size(x));
                for i=1:length(traitIdx)
                    s = obj.S*obj.scalars(ing);
                    t = traitIdx(i);
                    mean_posterior(:,i) = (1/( 1/s(t,t) + 1./sigmaEps(t,t) )).*(x(:,i)./sigmaEps(t,t));
                end
            end
            
        end
        
        % Calculate the p(alpha_hat|G=g)*p(G=g)
        function p_alpha_hat = p_alpha_hat(obj,x,sigmaEps,ing,traitIdx,marg)
            
            if ~marg
                p_alpha_hat = mvnpdf(x,zeros(size(x)),obj.S*obj.scalars(ing)+sigmaEps)*obj.ww(ing);
            else
                p_alpha_hat = zeros(size(x));
                for i=1:length(traitIdx)
                    t = traitIdx(i);
                    s = obj.S*obj.scalars(ing)+sigmaEps;
                    p_alpha_hat(:,i) = normpdf(x(:,i),0,sqrt(s(t,t)))*obj.ww(ing);
                end
            end
            
        end
        
        function a = get.F(obj)
            [a, ~] = obj.cptfn(obj.theta);
        end
        
        function a = get.dF(obj)
            [~, a] = obj.cptfn(obj.theta);
        end
        
        function a = get.noParams(obj)
            a = numel(obj.theta);
        end
        
        function a = get.noCons(obj)
            a = size(obj.vCon,2);
        end
        
        function c = get.cov(obj)
            c = obj.S * (obj.scalars * obj.ww);
        end
        
        function s = get.S(obj)
            [i,j]=find(triu(ones(obj.noTraits)));
            x2=zeros(obj.noTraits,length(i));
            x2(i' + (0:length(i)-1)*obj.noTraits)=1;
            x2(j' + (0:length(i)-1)*obj.noTraits) = ...
                x2(j' + (0:length(i)-1)*obj.noTraits)+1;
            x1=eye(obj.noTraits);
            b=obj.F(x1);
            s = 1/2 * (obj.F(x2) - b(i) - b(j));
            s = triuind(s);
            s = s + s' - diag(diag(s)); % symmetricize
        end
        
        function a = get.A(obj)
            if ~isempty(obj.vCon)
                a = obj.scalars.*obj.F(obj.vCon);
            else
                a = [];
            end
        end
        
        function a = get.dA(obj)
            if ~isempty(obj.vCon)
                a = (obj.scalars*obj.ww)*obj.dF(obj.vCon);
            else
                a = [];
            end
        end
        
        function a = get.noIngs(obj)
            a = numel([obj.scalars]);
        end
        
    end
end
