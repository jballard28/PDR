classdef ECF < matlab.mixin.Copyable
    properties
        traits {iscell} % trait names or file names
        T {ismatrix} % sampling times matrix
        phiBlocks {ismatrix} % ECF evaluated at sampling times within blocks
        phiEpsBlocks {ismatrix} % CF of sampling noise within blocks
        sigmaEps {ismatrix} % covariance matrix of sampling noise
        P {ismatrix} % projection/whitening matrix
    end
    
    properties (Dependent)
        phi % mean of phiBlocks
        phiPBlocks % P * phiBlocks
        phiP % P * phi
        noT % number of sampling times
        noTraits % number of traits
    end
    
    methods
        function obj = ECF(data,samplingtimes,varargin)
        % Class constructor
        %
        % Required Inputs:
        %   data: DATA object
        %   samplingtimes: vector of chosen sampling times
        %
        % Optional Inputs:
        %   whichTraits: Which traits from data to include. Can be vector
        %       of indices, or cell array of trait names. Order matters.
        %   maxNoSamplingtimes: Maximum number of sampling times
        %   noBlocks: Number of jackknife blocks
        %   tol: Tolerance for projection matrix. If this value is too
        %       small it may produce miscalibrated p-values.
        %
        % Outputs:
        %   obj: ECF object
        
            if nargin==0
                return;
            end
            
            % compute ecf of data and rotation matrix
            [obj.phiBlocks, obj.P, obj.sigmaEps, obj.T]=ecfMake(data,samplingtimes,varargin{:});

            %
            if ~isempty(data.traits)
                obj.traits = data.traits;
            else
                obj.traits = num2cell(1:obj.noTraits);
            end
        end
        
        function a = get.phi(obj)
        % Function to calculate and return the mean phi value over
        %   jackknife blocks
        %
        % Outputs:
        %   a: Mean phi over jackknife blocks
        
            a = mean(obj.phiBlocks,2);
            
        end
        
        function a = get.phiPBlocks(obj)
        % Function to get phi for each jackknife block after applying the
        %   projection matrix to prune correlated sampling times
        %
        % Outputs:
        %   a: phi for each jackknife block after pruning correlated
        %       sampling times
        
            a = obj.P*obj.phiBlocks;
            
        end
        
        function a = get.phiP(obj)
        % Function to get mean phi over jackknife blocks after applying the
        %   projection matrix to prune correlated samping times
        %
        % Outputs:
        %   a: Mean phi over jackknife blocks after pruning correlated
        %       sampling times
        
            a = mean(obj.phiPBlocks,2);
            
        end
        
        function a = get.noT(obj)
        % Function to get number of sampling time combinations
        %
        % Outputs:
        %   a: Number of sampling time combinations (before pruning) based
        %       on the number of ways to choose a sampling time for each trait
        %       from a vector of possible sampling times
        
            a = size(obj.T,2);
            
        end
        
        function a = get.noTraits(obj)
        % Function to calculate the number of traits
        %
        % Outputs:
        %   a: Number of traits
        
            a = size(obj.T,1);
            
        end
        
    end
    
    
end