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
        % class constructor
        function obj = ECF(data,samplingtimes,varargin)
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
            a = mean(obj.phiBlocks,2);
        end
        function a = get.phiPBlocks(obj)
            a = obj.P*obj.phiBlocks;
        end
        function a = get.phiP(obj)
            a = mean(obj.phiPBlocks,2);
        end
        function a = get.noT(obj)
            a = size(obj.T,2);
        end
        function a = get.noTraits(obj)
            a = size(obj.T,1);
        end
    end
    
    
end