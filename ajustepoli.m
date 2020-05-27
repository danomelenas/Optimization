clc
U = @(T) 204165.5./(330-2*T) + 10400./(T-20);
dx= @(x,delta) (U(x+delta)-U(x))./delta;
%ddx= @(x,delta) (U(x+delta)-2.*U(x)+ U(x-delta))./delta.^2;
delta=1e-4;
epsilon=1e-4;
x1=90;
x2=40;
fprintf('iter        x        F(x)      Fd(x) \n')
fprintf('----------------------\n')
tic
for i=1:100
    z=(3*(U(x1)-U(x2)))/(x2-x1)+dx(x1,delta)+dx(x2,delta);
    w=(x2-x1)/abs(x2-x1)*sqrt(z^2-dx(x1,delta)*dx(x2,delta));
    m=(dx(x2,delta)+w-z)/(dx(x2,delta)-dx(x1,delta)+2*w);
    if m<0
        xm=m;
    elseif m<=1
       xm=x2-m*(x2-x1);
    else
        xm=x1;
    end
   % fprintf('%d     %5.7f      %5.7f     %5.7f     \n',i,xm,U(xm),dx(xm,delta))
    if abs(dx(xm,delta))<epsilon
        break
    end
    if dx(xm,delta)*dx(x1,delta)<0
        x2=xm;
    else
        x1=xm;
    end
end
costo_c = toc;

fprintf('-------------------------\n')
fprintf('x* = %7.3f Mínimo = %8.3f\n',x,minimo)
fprintf('dx = %7.7f \n',abs(Dx(x,delta)))
fprintf('Número de iteraciones = %3d\n',i)
fprintf('Tiempo de procesamiento = %7.7f  (ms)\n',costo_c*1000)
