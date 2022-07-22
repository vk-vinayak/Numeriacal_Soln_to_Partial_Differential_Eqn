%This code solves one dimensional transient diffusion equation 
%( Parabolic partial differential equation) using theta method. It takes
%theta and time step (delta_t) as input. It gives temperature values at
%different space and time. Space X is discretized with delta_x=0.1 &
%temperature is calculated till 30 sec.This code also outputs a graph
%comparing variation in numeric error at different valus of X.
function a4q3
clc;
clear;
fprintf('\n Enter theta ');
theta=input('');
fprintf('\n Enter delta_t ');
delta_t=input('');
D=delta_t;
N=round(30/D);
%Defining initial and boundary conditions
T(1,:)=initial_f;
T(:,1)=0;
T(:,end)=0;
%Getting analytical soln 
y=zeros(N,11);
for j=1:11
y(1,j)=f(0,0.1*(j-1));
end
y(:,1)=0;
y(:,end)=0;
%Getting TDMA coefficients
for n=1:N-1
for i=1:9
    if (i==1)  
        b(i)=1+2*theta*D;
        c(i)=-theta*D;
        d(i)=(1-theta)*D*T(n,i+2)+(1-theta)*D*T(n,i)+((2*theta-2)*D+1)*T(n,i+1);
        
    elseif i==9 
        a(i)=-theta*D;
        b(i)=1+2*theta*D;
        d(i)=(1-theta)*D*T(n,i+2)+(1-theta)*D*T(n,i)+((2*theta-2)*D+1)*T(n,i+1);
    else
        a(i)=-theta*D;
        b(i)=1+2*theta*D;
        c(i)=-theta*D;
        d(i)=(1-theta)*D*T(n,i+2)+(1-theta)*D*T(n,i)+((2*theta-2)*D+1)*T(n,i+1);
    end
      y(n+1,i+1)=f((delta_t)*n,0.1*i) ;
         
end
  [b,d]=tdma(a,b,c,d);
  [x]=back(b,c,d);
  T(n+1,:)=[0 x 0];
  t_iter(n)=n;
end
T
err=T-y;
t_iter(N)=N;
%plotting Error graphs
plot(t_iter,err(:,6),'r-x',t_iter,err(:,4),'k-x',t_iter,err(:,2),'b-x');
title('Variation of Error with iteration');
legend('x=0.5','x=0.3','x=0.1');
xlabel('Time Iterations');
ylabel('Numerical error');
end
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
%Function to define initial condition
function T=initial_f
for i=1:11
 x=0.1*(i-1);
 T(i)=sin(pi*x);   
end
end
%Analytical Soln
function y=f(t,x)
alpha=0.01;
y=exp(-alpha*pi*pi*t)*sin(pi*x);
end
