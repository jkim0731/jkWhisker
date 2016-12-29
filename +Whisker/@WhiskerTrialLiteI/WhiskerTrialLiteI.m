classdef WhiskerTrialLiteI < Whisker.WhiskerTrialLite
    %
    %
    %
    %
    properties
        
        M0I = {};
        contactInds = {};
        
        

    end
    
    properties (Dependent = true)
        M0Combined
    end
    
    methods (Access = public)
        function obj = WhiskerTrialLiteI(w, varargin)
            %
            % USAGE:
            %
            %   obj = WhiskerSignalTrialI(w)
            %   obj = WhiskerSignalTrialI(w, 'polyRoiInPix',[roiMin, roiMax])
            %   obj = WhiskerSignalTrialI(w, 'polyRoiInPix',{[trajectoryIDs],[roiMin1, roiMax1], ...
            %                                     [roiMin2, roiMax2], ...
            %                                     [roiMin3, roiMax3], ...
            %                                     [roiMinN, roiMaxN]})
            %
            % INPUTS:
            %   w: a WhiskerTrial object.
            %
            %   Optional argument polyRoiInPix:
            %       Sets arc-length limits (in pixels) on which to perform secondary curve fitting.
            %       This argument can be given in two forms:
            %           (1) an 1x2 vector that gives the ROI for *all* whiskers; or
            %           (2) a cell array where first element is a vector of trajectory IDs
            %               (of length N) and subsequent elements comprise N 1x2 vectors
            %               giving ROIs for the trajectory IDs specified in the first
            %               element (respectively).
            %       Limits are inclusive.
            %   
            %   If WhiskerSignalTrial(w) is called without 'polyRoiInPix'
            %   argument, and w.polyFitsROI is empty (i.e.,
            %   fit_polys_roi() method of WhiskerTrial w was not called earlier), then
            %   an error is given.  A WhiskerSignalTrial requires
            %   polynomials to be fitted to the whiskers.
            %
            %
            
                      
            obj = obj@Whisker.WhiskerTrialLite(w,varargin{:});
           %   obj = obj@Whisker.WhiskerSignalTrial(w,varargin{:});
            
            if nargin==0
                return
            end
            
            obj.M0I = w.M0I;
            
            
        end
        
        
        
    end
    
    methods % Dependent property methods; cannot have attributes.
        
        function value = get.M0Combined(obj)
            
            if isempty(obj.barPos)
                value = [];
            else
               value = [1];
            end
        end
        
    end
    
end
