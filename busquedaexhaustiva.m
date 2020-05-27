%%%%%  exhaustive search or Brute-Force search
T=40:0.001:90;
U = @(T) 204165.5/(330-2*T) + 10400/(T-20);
delta=1e-4;
Dx= @(x,delta) (U(x+delta)- U(x))/ delta;
minimo=1e156;
x=0;
tic
for i=1:length(T)
    temp=U(T(i));
    if temp<minimo
        minimo=temp;
        x=T(i);
    end
end
costo_c = toc;
fprintf('-------------------------\n')
fprintf('x* = %7.3f Mínimo = %8.3f\n',x,minimo)
fprintf('dx = %7.7f \n',abs(Dx(x,delta)))
fprintf('Número de iteraciones = %3d\n',i)
fprintf('Tiempo de procesamiento = %7.7f  (ms)\n',costo_c*1000)