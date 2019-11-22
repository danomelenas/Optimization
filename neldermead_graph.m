function nlerdemead_graph

%%%%función objetivo
U = @(x)  100*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 +   90*(sqrt(x(1).^2  + (x(2)+1).^2) -1).^2 - (20.*x(1) + 40.*x(2));
U2 = @(x1,x2)  100*(sqrt(x1.^2  + (x2+1).^2) -1).^2 +   90*(sqrt(x1.^2  + (x2+1).^2) -1).^2 - (20.*x1 + 40.*x2);
[x1,x2]=meshgrid(-1:0.1:1,-1:0.1:1);


%%%% Condiciones iniciales
eps=1e-10; % condicion de paro 1D
xi = [-1 1];
limx1=[-1 1];
limx2=[-1 1];
x=xi;
fx_prev=U(x);

% rng('shuffle');
% x=rand(2,3);
% x(1,:)=x(1,:).*0.5*(limx1(2)-limx1(1))*(-1)^(randi([2,3],1)) ;
% x(2,:)=x(2,:).*0.5*(limx2(2)-limx2(1))* (-1)^(randi([2,3],1));
x_prev= [0.7143   0.1162;  -0.0292   0.365;   0.4006    0.8315]';
%
% x(2,:)=limx2(2)-x(2,:).*0.25*(limx2(2)-limx2(1));

color=[1,0,0;1,1,0;1,1,0];
fprintf('Initial function value = %7.4f\n ',fx_prev)
fprintf(' No. f(x) \n')
fprintf('__________________________________________\n')

for i=1:20
    %[x_curr,value]=graph_ite(x_prev,U,U2,color,x1,x2);
    [x_curr,value]=process_ite(x_prev,U);
    x_prev=x_curr;
    [~,order]=sort(value);
       
    fprintf('%3d %8.3f\n',i, U([x_curr(1,order(1)),x_curr(2,order(1))]))
end
end






function [x,value]=graph_ite(x,U,U2,color,x1,x2)
f=U2(x1,x2);
figure('color',[1 1 1]);
hold off
surf(x1,x2,f), shading interp
hold on
value=zeros(1,3);
punto=value;
alpha = 1;
gamma = 2;
rho=-0.5;
for i=1:3
    value(i)=U([x(1,i),x(2,i)]);
    punto(i)=plot3(x(1,i),x(2,i),value(i),'o','color',color(2,:), 'MarkerSize',10,'markerfacecolor',color(2,:));
    if i<3
        S=x(:,i+1)-x(:,i);
        for j= 0:0.01:1
            r=U(x(:,i)+j*S);
            xt=x(:,i)+j*S;
            plot3(xt(1),xt(2),r,'k.', 'MarkerSize',10,'markerfacecolor',[0 0 0])
        end
    else
        S=x(:,i)-x(:,1);
        for j= 0:0.01:1
            r=U(x(:,1)+j*S);
            xt=x(:,1)+j*S;
            plot3(xt(1),xt(2),r,'k.', 'MarkerSize',10,'markerfacecolor',[0 0 0])
        end
    end
    
end
[~,order]=sort(value);
posworst=order(end);



set(punto(posworst),'color',color(1,:))
set(punto(posworst),'markerfacecolor',color(1,:))
text(x(1,posworst),x(2,posworst),value(posworst)+20, 'x_{worst}','FontSize',18)

% Centroid
xc=(sum(x,2)-x(:,posworst))/2;
plot3(xc(1),xc(2),U([xc(1),xc(2)]),'o','color',[1 0 1], 'MarkerSize',10,'markerfacecolor',[1 0 1])
text(xc(1),xc(2),U([xc(1),xc(2)])+20, 'x_{c}','FontSize',18)

% Reflection
xr=xc+alpha *(xc-x(:,posworst));
plot3(xr(1),xr(2),U([xr(1),xr(2)]),'o','color',[0.5 1 0.5], 'MarkerSize',10,'markerfacecolor',[0.5 1 0.5])
text(xr(1),xr(2),U([xr(1),xr(2)])+30, 'x_{r}','FontSize',18)



if U([xr(1),xr(2)])< value(order(end-1)) &&  U([xr(1),xr(2)])> value(order(1))
    x(:,posworst)=xr;
elseif  U([xr(1),xr(2)]) < value(order(1))
    
    % Expansion
    xe = xc + gamma*(xc-x(:,posworst));
    plot3(xe(1),xe(2),U([xe(1),xe(2)]),'o','color',[0.2 0.5 0.7], 'MarkerSize',10,'markerfacecolor',[0.2 0.5 0.7])
    text(xe(1),xe(2),U([xe(1),xe(2)])+30, 'x_{e}','FontSize',18)
    
    if U([xe(1),xe(2)])< U([xr(1),xr(2)])
        x(:,posworst)=xe;
    else
        x(:,posworst)=xr;
    end
elseif U([xr(1),xr(2)]) > value(order(end-1))
    % Contraction
    xcon = xc + rho*(xc-x(:,posworst));
    plot3(xcon (1),xcon (2),U([xcon(1),xcon(2)]),'o','color',[0.2 0.5 0.7], 'MarkerSize',10,'markerfacecolor',[0.2 0.5 0.7])
    text(xcon (1),xcon (2),U([xcon(1),xcon(2)])+30, 'x_{e}','FontSize',18)
    
    if U([xcon(1),xcon(2)])<  value(order(end))
        x(:,posworst)=xcon;
    else
        for i=2:3
            x(:,order(i))=x(:,1)-rho*(x(:,order(i))-x(:,1));
        end
    end
    
end

end






function [x,value]=process_ite(x,U)
value=zeros(1,3);
alpha = 1;
gamma = 2;
rho=-0.5;
for i=1:3
    value(i)=U([x(1,i),x(2,i)]);
end
[~,order]=sort(value);
posworst=order(end);

% Centroid
xc=(sum(x,2)-x(:,posworst))/2;

% Reflection
xr=xc+alpha *(xc-x(:,posworst));


if U([xr(1),xr(2)])< value(order(end-1)) &&  U([xr(1),xr(2)])> value(order(1))
    x(:,posworst)=xr;
elseif  U([xr(1),xr(2)]) < value(order(1))
    % Expansion
    xe = xc + gamma*(xc-x(:,posworst));
    
    if U([xe(1),xe(2)])< U([xr(1),xr(2)])
        x(:,posworst)=xe;
    else
        x(:,posworst)=xr;
    end
elseif U([xr(1),xr(2)]) > value(order(end-1))
    % Contraction
    xcon = xc + rho*(xc-x(:,posworst));
    if U([xcon(1),xcon(2)])<  value(order(end))
        x(:,posworst)=xcon;
    else
        for i=2:3
            x(:,order(i))=x(:,1)-rho*(x(:,order(i))-x(:,1));
        end
    end
    
end

end