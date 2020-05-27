U = @(T) 204165.5./(330-2*T) + 10400./(T-20);
dx= @(x,delta) (U(x+delta)-U(x))./delta;
ddx= @(x,delta) (U(x+delta)-2.*U(x)+ U(x-delta))./delta.^2;
 delta=1e-4;
 epsilon=1e-4;
xc=90;
fprintf(' iter  x        F(x)             Fd(x)          \n ')

fprintf('-----------------------------------------\n ')
tic
for i=1:100
  
    xp=xc;
    xc=xp - dx(xp, delta)/ddx(xp,delta);
    if abs (dx(xc,delta)) < epsilon
        break
    end
   
end
 costo_c = toc;
   
 fprintf('-----------------------------------------\n ')

   fprintf('%d        %5.7f      %5.7f     %5.7f   \n ',i,xc,U(xc),dx(xc,delta));
  fprintf('Número de iteraciones = %3d\n',i)
fprintf('Tiempo de procesamiento = %7.7f  (ms)\n',costo_c*1000)
