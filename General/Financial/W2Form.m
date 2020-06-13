classdef W2Form
    % This class stores the information contained on a W2 tax form.
    
    properties
        Boxes = TaxFormField(Name, Value, Description);
    end
    
    methods
        function obj = W2Form()
            
            box_names = {'a', 'b', 'c', 'd', 'e', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12a', '12b', '12c', '12d', '13', '14', '15', '16', '17', '18', '19', '20'};
            box_values = {'', '', '', [], ''};
            
        end
        
%         function obj = untitled6(inputArg1,inputArg2)
%             %UNTITLED6 Construct an instance of this class
%             %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
%         end
%         
%         function outputArg = method1(obj,inputArg)
%             %METHOD1 Summary of this method goes here
%             %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
%         end
    end
end

