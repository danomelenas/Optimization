function newton_multivariable

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);



%%%% Condiciones iniciales
delta=1e-3; % para derivada
ep1=1e-5; % condicion de paro 1D
ep2=1e-5; % condicion de paro 2D
xi = [-1 1];
x = xi;
fprintf('Initial function value = %7.4f\n',U(x))
fprintf('No. x-vector f(x) Deriv \n')
fprintf('__________________________________________\n')

for i = 1:50
fx_prev=U(x);
grad(1) = Dx1(x,delta);
grad(2) = Dx2(x,delta);
H = matrizH(x,delta,U);

x = (x' - H\grad')';
f = U(x);
if abs(f-fx_prev)<ep1 || norm(grad)<ep2
break;
end
fprintf('%3d %8.3f %8.3f % 8.3f %8.3f\n',i,x,f,norm(grad))
end
fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,f,norm(grad))
fprintf('__________________________________________\n')

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