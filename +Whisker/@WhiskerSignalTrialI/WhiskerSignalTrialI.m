classdef WhiskerSignalTrialI < Whisker.WhiskerSignalTrial
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
        function obj = WhiskerSignalTrialI(w, varargin)
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
            
            extrapDistanceInPix = 32;
            roi = [64 3*32];
                      
            obj = obj@Whisker.WhiskerSignalTrial(w,varargin{:});
            
            if nargin==0
                return
            end
            
            obj.recompute_cached_follicle_coords(extrapDistanceInPix);
            
            for k=1:length(obj.trajectoryIDs)
                obj.M0I{k} = obj.calc_M0I(obj.trajectoryIDs(k),roi);
            end
            
        end
      

        
        function r = calc_M0I(obj,tid,roi)
           disp('Calculating inertial force.') 
           [th,t] = obj.get_fitted_line_angle_in_roi(tid,roi);
           
            m=18.8e-9 % Mass of whisker in kilograms
            h=16e-3 % length of whisker in m
            r=33.5e-6 % radius of whisker at base in m
            
           I=3/20*m*(r^2+4*h^2)+m*h^2/16;
           r = I*([0 0 diff(diff(smooth(deg2rad(th),3)))']./[0 diff(obj.get_time(tid))].^2);
           
        end
        
        
        
    end
    
    methods % Dependent property methods; cannot have attributes.
        
        function value = get.M0Combined(obj)
            
            if isempty(obj.barPos)
                value = [];
            else
               value = [];
            end
        end
        
    end
    
end
