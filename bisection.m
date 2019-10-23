function bisection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB code bisection.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% delx -> x for central difference method

clear
clc
U = @(T) 204165.5/(330-2*T) + 10400/(T-20);
a = 40; %-> lower bound of the design variable
b = 90; %-> upper bound of the design variable
delx = 0.01; %-> x for central difference method
epsilon = 0.001; %->  stop threshold
fprintf('     a       b           f(x)\n')
fprintf('----------------------------------\n')
for i= 1:100
    fprintf(' %7.3f %8.3f      %8.3f\n',a,b,U(a))
    alpha = (a+b)/2;
    derivative = (U(a+delx) - U(a-delx) )/(2*delx);
    derivative_alpha = (U(alpha+delx)- U(alpha-delx))/(2*delx);
    if (derivative*derivative_alpha) < 0
        b = alpha;
    else
        a = alpha;
    end
    if abs(a-b) < epsilon
        break;
    end
end
fprintf('----------------------------------\n')
fprintf('x* = %7.3f Minimum = %8.3f\n',a,U(a))
fprintf('Number of function calls = %3d\n',i)
end

