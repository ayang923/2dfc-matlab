classdef curve_param_obj < handle
    %CURVE_PARAM_OBJ Struct like object for storing discretization
    %parameters for curve
    
    properties
        curve_n
        n_total
        start_idx
        end_idx
        
        U_c_0_intervals
        U_c_1_intervals
    end
    
    methods
        function obj = curve_param_obj(curve_n)
            %CURVE_PARAM_OBJ Construct an instance of this class
            %   Detailed explanation goes here
            obj.curve_n = curve_n;
            obj.n_total = sum(curve_n);
            obj.start_idx = cumsum([1, curve_n(1:end-1)'])';
            obj.end_idx = obj.start_idx+curve_n-1;
            
            obj.U_c_0_intervals = zeros(length(curve_n), 2);
            obj.U_c_1_intervals = zeros(length(curve_n), 2);
        end
    end
end

