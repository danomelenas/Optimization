function levenbermarquardt

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);



%%%% Condiciones iniciales
delx=1e-3; % para derivada
eps=1e-5; % condicion de paro 1D
xi = [-1 1];
lambda = 1e3;
x = xi;
fprintf('Initial function value = %7.4f\n',U(x))
fprintf('No. x-vector f(x) Deriv \n')
fprintf('__________________________________________\n')
for i = 1:100
    fx_prev=U(x);
    grad(1) = Dx1(x,delx);
    grad(2) = Dx2(x,delx);
    H = matrizH(x,delx,U);
    search = -inv(H+lambda*eye(length(x)))*grad';
    x = x + search';
    f = U(x);
    if f < fx_prev
        lambda = lambda/2;
    else
        lambda = 2*lambda;
    end
    if abs(f-fx_prev)<eps || norm(grad)<eps
        break;
    end
    fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,f,norm(grad))
end
fprintf('__________________________________________\n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



function [H]=matrizH(x,delta,U)
dt=[delta,0];
H(1,1)= (U(x+dt)- 2*U(x) + U(x-dt))/ dt(1).^2;
dt=[0,delta];
H(2,2)= (U(x+dt)- 2*U(x) + U(x-dt))/ dt(2).^2;
H(2,1)= (U([x(1)+delta,x(2)+delta])  -  U([x(1)+delta,x(2)-delta]) -...
    U([x(1)-delta,x(2)+delta])  +  U([x(1)-delta,x(2)-delta]))/...
    (4*delta^2);
H(1,2)= H(2,1);

end