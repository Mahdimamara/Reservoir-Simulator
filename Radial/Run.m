% Run the Main program
clear
clc
fprintf('1:Linear flow regime\n');
fprintf('2:Radial flow regime\n');
t_fl=input('Type of flow regime: ');
switch t_fl
    case{1}
        run linear_RM
    case{2}
        run Radial_RM
end