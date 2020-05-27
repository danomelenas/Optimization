U = @(T) 204165.5./(330-2*T) + 10400./(T-20);
dx= @(x,delta) (U(x+delta)-U(x))./delta;
%ddx= @(x,delta) (U(x+delta)-2.*U(x)+ U(x-delta))./delta.^2;
delta=0.01;
epsilon=0.01;
b=90;
a=40;
fprintf('iter        x        F(x)      Fd(x) \n')
fprintf('----------------------\n')
for i=1:100
   lw=b-a;
   w1=a+0.6180339887*lw;
   w2=b-0.6180339887*lw;
    fx1=U(w1);
    fx2=U(w2);
 if fx1<fx2
    a=w2;
 else
     b=w1;
 end
    lw=b-a;
    r=(a+b)/2;
    fprintf('%d     %5.7f      %5.7f     %5.7f     \n',i,r,U(r),dx(r,delta))

if abs(lw)<epsilon
    break
end
end