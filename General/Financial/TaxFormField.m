classdef TaxFormField
    % This class stores the information contained in a tax form field.
    
    properties
        Name = '';
        Value = [];
        Description = '';
    end
    
    methods
        
        % Implement this object's constructor.
        function obj = TaxFormField(Name, Value, Description)
            
            % Assign this objects property values.
            if nargin >= 1, obj.Name = Name; end
            if nargin >= 2, obj.Value = Value; end
            if nargin >= 3, obj.Description = Description; end

        end
        
    end
end

