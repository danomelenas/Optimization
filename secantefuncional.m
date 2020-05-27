clc
U = @(T) 204165.5./(330-2*T) + 10400./(T-20);
dx= @(x,delta) (U(x+delta)-U(x))./delta;
delta=1e-4;
epsilon=1e-4;
x1=90;
x2=40;
%fprintf('iter     x     F(x)      Fd(x) \n')
fprintf('----------------------\n')
tic
for i=1:100
    alpha=x2-(dx(x2,delta)*(x2-x1))/(dx(x2,delta)-dx(x1,delta));
    if dx(alpha,delta)*dx(x1,delta)<0
        x2=alpha;
    else
        x1=alpha;
    end
    
    if abs(dx(alpha,delta))<epsilon
        break
    end
end
fprintf('%d     %5.7f      %5.7f     %5.7f     \n',i,alpha,U(alpha),dx(alpha,delta))   
 costo_c = toc;

fprintf('-------------------------\n')

fprintf('Tiempo de procesamiento = %7.7f  (ms)\n',costo_c*1000)
