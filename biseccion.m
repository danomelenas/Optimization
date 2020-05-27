Uf=@(T) 204165.5./(330-2.*T)+10400./(T-20);
Dx=@(x,delta) (Uf(x+delta)-Uf(x))/delta;
a=40;
b=90;
epsilon=1e-4;
delta=0.01;
tic
for i=1:5000
    alpha=a+(b-a)/2;
    fpx=Dx(alpha,delta);
if Dx(a,delta)*fpx<0
     b=alpha;
else
      a=alpha;
end

if  abs(fpx)<epsilon
  
    break
end

end

costo_c = toc;
fprintf('valor final de la funcion U= %7.5f T=%2.3f  deri=%5.5f \n',Uf(alpha),alpha,fpx)


fprintf('-------------------------\n')
fprintf('x* = %7.3f Mínimo = %8.3f\n',x,minimo)
fprintf('dx = %7.7f \n',abs(Dx(x,delta)))
fprintf('Número de iteraciones = %3d\n',i)
fprintf('Tiempo de procesamiento = %7.7f  (ms)\n',costo_c*1000)
