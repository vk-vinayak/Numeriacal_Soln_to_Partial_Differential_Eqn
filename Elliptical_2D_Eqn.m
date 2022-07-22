%This program solves two dimensional elliptical partial differential equation in T
%(Temperature). Physics-wise,it is 2-D thermal conduction equation with temperature boundary
%condition for the rectangular slab which is 1m in x direction & 2m in y direction. 
% The code uses point SOR method with w as over relaxation parameter.
%Over relaxation parameter w can be varied between 1 to 2 to minimize
%number of iteration. For this grid size of 21 x 41, optimum w as 1.785.
%This code outputs matrices outputting analytical (matrix Y) and
%numerical (matrix T) solution. It also outputs 2 graphs which compare
%analytical & numerical solutions and numerical error associated with this
%method.
function a5b
clc;
clear;
L=1;
H=2;
n=21;
m=41;
x=linspace(0,L,n);
y=linspace(0,H,m);
delta_x=L/(n-1);
delta_y=H/(m-1);
b=delta_x/delta_y;
w=1.785;
T=zeros(m,n);
T(1,:)=0;
T(:,1)=0;
T(:,end)=0;
T(end,:)=1;
M=zeros(m,n);
M(1,:)=0;
M(:,1)=0;
M(:,end)=0;
M(end,:)=1;
%Analytical soln
for k=1:m
for l=1:n
x=delta_x*(l-1);
y=delta_y*(k-1);
Y(k,l)=analytical(x,y,L,H);
end
end
%Numeric soln
add=1;
v=0;
while add >0.01
for j=m-1:-1:2
for i=n-1:-1:2
M(j,i)=(1-w)*T(j,i)+w*(T(j,i+1)+T(j,i-1)+(b^2)*(T(j-1,i)+T(j+1,i)))/(2*(1+b^2));
err(j,i)=M(j,i)-T(j,i);
T(j,i)=M(j,i);
end
end
v=v+1;
add=sum(sum(abs(err)));
end
%Plotting graphs at x=0.5
u=linspace(0,H,m);
plot(u,T(:,(end+1)/2),'b-o',u,Y(:,(end+1)/2),'r-x');
title('Variation of Temperature at x=0.5 in Vertical Direction');
legend('Numerical','Analytical');
xlabel('y');
ylabel('Temperature');
figure
plot(u,T(:,(end+1)/2)-Y(:,(end+1)/2),'g-o');
title('Variation of Error(Numeric-Analytical) at x=0.5 in Vertical Direction');
xlabel('y');
ylabel('Error');
%T is numeric solution
T
%Y is analytical solution
Y

end
function term=analytical(x,y,L,H)
term=0;
for n=1:20
   term=term+((1+(-1)^(n+1)))*(sinh(n*pi*(y)/L))*(sin(n*pi*x/L))/(sinh(n*pi*H/L)*n*pi) ;
end
term=term*2;
end