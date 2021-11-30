classdef COMPONENT < matlab.mixin.Copyable
    properties
        theta % column vector of parameters for this component
        scalars % one scalar for each ingredient
        name % name of component (just for convenience)
        ww % weight of each ingredient
        thetaCon % upper and lower bounds on theta
        noTraits % number of traits
        grad % gradient of the objective function
        conType % constraint type ('none' or 'diagonal' - see h2ConOpt below)
        fixedParam % parameter (different from theta) that isn't learned
    end
    
    properties (Dependent)
        A % matrix of constraint values, no. constraints x no. ingredients
        dA % gradient of A
        noIngs % length(scalars)
        noCons % number of coefficient constraints
        noParams % number of parameters
        cov % S*dot(scalars,ww)
        S % covariance matrix of the component
    end
    
    properties (Access = private)
        cptfn % used to compute gradients
        vCon % used to compute constraint coefficients
        type % component type - used to determine which component functions to use
    end
    
    properties (Dependent, Access = private)
        F % component function
        dF % derivative of component function
    end
    
    methods
        function obj = COMPONENT(type,no_traits,scalars,fixedParam,varargin)
        % Class constructor
        %
        % Required Inputs:
        %   type: Which type of component (member of cpts folder)('full_rank','rank_one','fixed')
        %   no_traits: Number of traits
        %   scalars: Scalars, one per ingredient
        %   fixedParam: Parameter (different from theta) that isn't learned
        %
        % Optional Inputs:
        %   theta: Parameters; if unspecified, these are randomly initialized
        %   ww: Weight of each ingredient
        %   name: Name of component (just for convenience)
        %   h2ConOpt: option to generate heritability constraints. Options allowed:
        %       'none' (no constraint), 'diagonal' (constrain h2, but not
        %       genetic correlation)
        %
        % Outputs:
        %   obj: COMPONENT object
        
            if nargin==0
                return;
            end
            p=inputParser;
            
            addRequired(p, 'type', @isstr);
            addRequired(p, 'no_traits', @(x)isnumeric(x) & isscalar(x));
            addRequired(p, 'scalars', @(x)isnumeric(x) & isvector(x) & size(x,1)==1);
            addRequired(p, 'fixedParam', @(x)true)
            
            [~, ~, randomParams, obj.thetaCon] = component_function(type,no_traits,fixedParam,1);
            
            addOptional(p, 'theta', randomParams, ...
                @(x)isnumeric(x) & all(size(x)==size(randomParams)));
            addOptional(p, 'ww', zeros(size(scalars))', ...
                @(x)isnumeric(x) & size(x,1) == size(scalars,2));
            addOptional(p, 'name', type, @isstr)
            addOptional(p, 'h2ConOpt', 'none', @isstr)
            
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
            i=strcmp(p.Results.h2ConOpt,{'none','diagonal'});
            
            if i(1)
                obj.vCon=zeros(no_traits,0);
            elseif i(2)
                obj.vCon=eye(no_traits);
            else
                error('h2ConOpt must be either none or diagonal')
            end
            
            obj.vCon=zeros(no_traits,0);
            obj.scalars = scalars;
            obj.ww = p.Results.ww;
            obj.name = p.Results.name;
            
        end
        
        function newCon(obj,v)
        % Function to add a constraint
            obj.vCon = [obj.vCon;v];
        end
        
        function delCon(obj,kk)
        % Function to delete a constraint
        %
        % Optional Inputs:
        %   kk: Index of constraint to delete; otherwise deletes all
        %       constraints
        
            if nargin < 1
                kk = 1:obj.noCons;
            end
            obj.vCon = obj.vCon(setdiff(1:end,kk),:);
        end
        
        function settheta(obj,thetaNew)
        % Function to set theta to a new value while respecting parameter constraints
        %
        % Inputs:
        %   thetaNew: New theta values
        
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
        
        function dX = gradCompute(obj,ecf,xbeta,additive_cpts)
        % Function to evaluate the gradient of the RSS wrt theta
        %
        % Required Inputs:
        %   ecf: ECF object
        %   xbeta: Regression weights multiplied by characteristic function
        %       values
        %
        % Optional Inputs:
        %   additive_cpts: Whether additive components are being used
        %
        % Outputs:
        %   dX: Gradient of the objective function wrt theta
        
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
                
                % multiply gradient by projection matrix
                dX = ecf.P * dX;
                
                % gradient of the objective function
                obj(i).grad = -2 * (ecf.phiP - ecf.P * (xbetaProd - 1))'*dX;
            end
        end
        
        function randomize(obj)
        % Function to randomize the parameters of the component
        
            for i=1:numel(obj)
                [~, ~, randomParams] = component_function(obj(i).type,obj(i).noTraits,obj(i).fixedParam);
                obj(i).settheta(randomParams);
            end
        end
        
        function b=simulate(obj,mm)
        % Function to simulate from a component
        %
        % Inputs:
        %   mm: Number of SNPs
        %
        % Outputs:
        %   b: Vector of effect sizes simulated from the component
        
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
        
        function step(obj, stepsize)
        % Function to take a step down the gradient (supports nonscalar objects)
        %
        % Inputs:
        %   stepsize: Size of the step to take down the gradient
        
            for i = 1:numel(obj)
                obj(i).settheta(obj(i).theta - stepsize * obj(i).grad');
            end
        end
        
        function setww(obj, wwNew)
        % Function to set ww to a new value (supports nonscalar objects)
        %
        % Inputs:
        %   wwNew: Values to set as new weights
        
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
        % Function to calculate the characteristic function value at
        %   specified sampling times
        %
        % Inputs:
        %   input: Either an ECF object or vector of sampling times
        %
        % Outputs:
        %   a: Objective function value at each of the sampling times for
        %       each mixture Gaussian (corresponding to a given scalar)
        
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
        
        function a = get.F(obj)
        % Get the function for the component
        %
        % Outputs:
        %   a: Component function
        
            [a, ~] = obj.cptfn(obj.theta);
            
        end
        
        function a = get.dF(obj)
        % Function to get the derivative of the component function wrt
        %   theta
        %
        % Outputs:
        %   a: Derivative of the component function wrt theta
        
            [~, a] = obj.cptfn(obj.theta);
            
        end
        
        function a = get.noParams(obj)
        % Function to get the number of parameters in a component
        %
        % Outputs:
        %   a: Number of parameters
        
            a = numel(obj.theta);
            
        end
        
        function a = get.noCons(obj)
        % Function to get the number of constraints in the component
        %
        % Outputs:
        %   a: Number of constraints
        
            a = size(obj.vCon,2);
            
        end
        
        function c = get.cov(obj)
        % Function to get the covariance matrix for a component
        %
        % Outputs:
        %   cov: Component covariance matrix (pattern matrix multiplied by
        %       expected scaling parameter)
        
            c = obj.S * (obj.scalars * obj.ww);
            
        end
        
        function s = get.S(obj)
        % Function to get the pattern matrix for the component
        %
        % Outputs:
        %   s: Component pattern matrix
        
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
        % Function to get the matrix of constraint values
        %
        % Outputs:
        %   a: Constraint values
        
            if ~isempty(obj.vCon)
                a = obj.scalars.*obj.F(obj.vCon);
            else
                a = [];
            end
            
        end
        
        function a = get.dA(obj)
        % Function to get the gradient of A
        %
        % Outputs:
        %   a: Gradient of A
        
            if ~isempty(obj.vCon)
                a = (obj.scalars*obj.ww)*obj.dF(obj.vCon);
            else
                a = [];
            end
            
        end
        
        function a = get.noIngs(obj)
        % Function to get the number of ingredients (number of mixture
        %   Gaussians within the component)
        %
        % Outputs:
        %   a: Number of ingredients in the component
        
            a = numel([obj.scalars]);
            
        end
        
    end
end
