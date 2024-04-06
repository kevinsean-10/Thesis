classdef MyClass < handle
    properties
        Property1
        Property2
    end
    
    methods
        function obj = MyClass(value1, value2)
            % Constructor method to initialize properties
            obj.Property1 = value1;
            obj.Property2 = value2;
        end
        
        function updateProperties(obj, newValue1, newValue2)
            % Method to update properties
            obj.Property1 = newValue1;
            obj.Property2 = newValue2;
        end
    end
end
