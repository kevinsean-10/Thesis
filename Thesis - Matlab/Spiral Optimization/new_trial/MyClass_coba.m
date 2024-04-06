clear;clc

% Create an instance of myClass
object = MyClass(2, 3);

% Display the initial properties
disp(object.Property1);  % Output: 2
disp(object.Property2);  % Output: 3

% Call the updateProperties method to modify the properties
object.updateProperties(4, 6);

% Display the updated properties
disp(object.Property1);  % Output: 4
disp(object.Property2);  % Output: 6

dodo = @(x,y) x + y(1);
dodo(1,[2,3]);

% [[1,2],[-1,5]]
% % x1>1 % x1-1>0
% % x2>-1 % x2-(-1)>0
% % -x1>-2 % -(x1-2)>0
% % -x2>-5 % -(x2-5)>0

