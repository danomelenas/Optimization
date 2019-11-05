function CG

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));

%%%% Derivadas parciales para gradiente
Dx1= @(x,delta)   (U([x(1)+delta,x(2)])- U([x(1)-delta,x(2)]))/ (2*delta);
Dx2= @(x,delta)   (U([x(1),x(2)+delta])- U([x(1),x(2)-delta]))/ (2*delta);



%%%% Condiciones iniciales
delx=1e-3; % para derivada
eps=1e-3; % condicion de paro 1D
xi = [-1 1];
x = xi;
fprintf('Initial function value = %7.4f\n',U(x))
fprintf('No. x-vector f(x) Deriv \n')
fprintf('__________________________________________\n')
for i = 1:300
    if i==1
        grad_prev(1) = Dx1(x,delx);
        grad_prev(2) = Dx2(x,delx);
        search_prev = -grad_prev;
        [alpha,fx] = seccion_dorada(x,search_prev ,xi,eps,U);
        
        
        if norm(grad_prev)<eps
            break;
        end
        x = x + alpha*search_prev;
        fx_prev = func_multivar(x);
    else
        grad(1) = Dx1(x,delx);
        grad(2) = Dx2(x,delx);
        search = -grad +((norm(grad)^2)/(norm(grad_prev)^2))*search_prev;
        [alpha,fx] = seccion_dorada(x,search,xi,eps,U);
        if abs(fx-fx_prev)<eps ||norm(grad)<eps
            break;
        end
        grad_prev = grad;
        search_prev = search;
        x = x + alpha*search;
        fx_prev = U(x);
    end
    fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,fx,norm(grad_prev))
end
fprintf('%3d %8.3f %8.3f % 8.3f %8.3f \n',i,x,fx,norm(grad))
fprintf('__________________________________________\n')
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
