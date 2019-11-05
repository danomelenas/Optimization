function modified_newton

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);
xi = [-1 1];
x = xi;
eps= 1e-3;
delx = 1e-3;
fx_prev=U(x);
fprintf('Initial function value = %7.4f\n ',fx_prev)
fprintf('No. x-vector f(x) Deriv \n')
fprintf('__________________________________________\n')
for i = 1:50
    grad(1) = Dx1(x,delx);
    grad(2) = Dx2(x,delx);
    H = matrizH(x,delx,U);
    search = -inv(H)*grad';
    [alpha,f ] = seccion_dorada(x,search',xi,eps,U);
    if abs(f -fx_prev)<eps || norm(grad)<eps
        break;
    end
    fx_prev = f;
    x = x + alpha*search';
    fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,f,norm(grad))
end
fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,f,norm(grad))
fprintf('__________________________________________\n')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [alpha1,falpha1]=seccion_dorada (x,Si,lims,ep1,U)

tau = 0.381967;
alpha1 = lims(1)*(1-tau) + lims(2)*tau;
alpha2 = lims(1)*tau + lims(2)*(1-tau);
falpha1 = U(x+alpha1*Si);
falpha2 = U(x+alpha2*Si);

for i= 1:1000
    if falpha1 > falpha2
        lims(1) = alpha1;
        alpha1 = alpha2;
        falpha1 = falpha2;
        alpha2 = tau*lims(1) + (1-tau)*lims(2);
        falpha2 = U(x+alpha2*Si);
    else
        lims(2) = alpha2;
        alpha2 = alpha1;
        falpha2 = falpha1;
        alpha1 = tau*lims(2) + (1-tau)*lims(1);
        falpha1 = U(x+alpha1*Si);
    end
    if abs(U(x+alpha1*Si)- U(x+alpha2*Si)) < ep1
        break;
    end
end
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