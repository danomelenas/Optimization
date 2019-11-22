function powell

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));


%%%% Condiciones iniciales
eps=1e-8; % condicion de paro 1D
xi = [-1 1];
x=xi;
fx_prev=U(x);
fprintf('Initial function value = %7.4f\n ',fx_prev)
fprintf(' No. x-vector f(x) \n')
fprintf('__________________________________________\n')


for i = 1:2
    for j = 1:3
        if (i==j)
            term(i,j)=1;
        else
            term(i,j) = 0;
        end
    end
end

for i = 1: 2
    search{i} = (term(:,i))';
end


for ite = 1:100
    xini = x;
    i = 1;
    while i<3
        [alpha,fx] = seccion_dorada (x,search{i},xi,eps,U);
        x = x + alpha*search{i};
        i = i+1;
    end
    
    if abs(fx-fx_prev) < eps
        break;
    end
    search{i} = (x-xini);
    [alpha,fx] = seccion_dorada (x,search{i},xi,eps,U);
    x = x + alpha*search{i};
    temp = search;
    for i = 1:2
        search{i} = temp{i+1};
    end
    fx_prev = fx;
    fprintf('%3d %8.3f %8.3f %8.3f\n',ite, x ,fx_prev)
    
end
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