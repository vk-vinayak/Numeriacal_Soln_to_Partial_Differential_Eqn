%This code solves two dimensional transient diffusion equation (Parabolic
%partial differential equation). It uses fractional step ADI method. It
%takes input like number of x & y grids (10), time (4sec) & delta_t (0.1). It outputs
%temperature matrix and 2 graphs for analytical (Figure 1) & numerical
%solution (Figure 2).
function a6q11
alpha=10^-2;
aa=input('Enter no. of x grid');
bb=input('Enter no. of y grid');
u=input('Enter time untill you want analysis');
x_grid=1/aa;
y_grid=1/bb;
t_grid=input('Enter delta t');
n=1/(x_grid)+1;
m=1/(y_grid)+1;
p=1/(t_grid);
bx=alpha.*t_grid/(2*x_grid.^2);
by=alpha.*t_grid/(2*y_grid.^2);
T=zeros(m,n);
e=1;
k=1;
%Initial Condition
for j=1:m
 for i=1:n
    x=(i-1).*x_grid;
    y=(j-1).*y_grid;
    T(j,i)=sin(pi*y).*sin(pi*x);
 end
end

for k=1:t_grid:u+1-t_grid
%Boundary Conditions
T(1,:)=0;
T(end,:)=0;
T(:,1)=0;
T(:,end)=0;
for j=1:m
    for i=1:n
        x=(i-1).*x_grid;
        y=(j-1).*y_grid;
        z=(k-1).*t_grid;
        analytic(j,i)=sin(pi*y).*sin(pi*x).*exp(-2*0.01*pi*pi*z);
    end
end
%Calculation from 0 to n+1/2 time step
for p=1:n-2
for i=1:n-2
    if (i==1)  
        b(i)=-(2*bx+1);
        c(i)=bx;
        d(i)=-by*(T(p,i+1)+T(p+2,i+1))+(2*by-1)*T(p+1,i+1);
        
    elseif i==n-2 
        a(i)=bx;
        b(i)=-(2*bx+1);
        d(i)=-by*(T(p,i+1)+T(p+2,i+1))+(2*by-1)*T(p+1,i+1);
    else
        a(i)=bx;
        b(i)=-(2*bx+1);
        c(i)=bx;
        d(i)=-by*(T(p,i+1)+T(p+2,i+1))+(2*by-1)*T(p+1,i+1);
    end  
end 
  [b,d]=tdma(a,b,c,d);
  [x]=back(b,c,d);
  T(p+1,:)=[0 x 0];
 %Calculation from n+1/2 to n+1 time step
for j=1:m-2
    if (j==1)  
        B(j)=-(2*by+1);
        C(j)=by;
        D(j)=-bx*(T(j+1,p)+T(j+1,p+2))+(2*bx-1)*T(j+1,p+1);
        
    elseif j==m-2 
        A(j)=by;
        B(j)=-(2*by+1);
        D(j)=-bx*(T(j+1,p)+T(j+1,p+2))+(2*bx-1)*T(j+1,p+1);
    else
        A(j)=by;
        B(j)=-(2*by+1);
        C(j)=by;
        D(j)=-bx*(T(j+1,p)+T(j+1,p+2))+(2*bx-1)*T(j+1,p+1);
    end
end 
  [B,D]=tdma(A,B,C,D);
  [y]=back(B,C,D);
  T(:,p+1)=[0 y 0];
end

end
 T
 
 
 figure
 contourf(0:x_grid:1 ,0:y_grid:1 ,T)
 figure
 contourf(0:x_grid:1 ,0:y_grid:1 ,analytic)
end
 %Function for thomas algorithm
function [b,d]=tdma(a,b,c,d)
n=length(b);
for i=2:n
    m=a(i)/b(i-1);
    b(i)=b(i)-m*c(i-1);
    d(i)=d(i)-m*d(i-1);
end
end
function [x]=back(b,c,d)
n=length(b);
x(n)=d(n)/b(n);
for i=n-1:-1:1
    x(i)=(d(i)-c(i)*x(i+1))/b(i);
end
end