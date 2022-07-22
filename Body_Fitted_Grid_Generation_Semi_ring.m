%This code generates body fitted laplacian grid (NXN)for semi ring with
%inner radius as 1 unit & outer radius as 5 units for 11 x 11 grid. Inner
%radius temperature is 30 while the outside is 100 units. The code outputs
%analytical and numerical temperature contour along with generated body
%fitted grids. It also outputs dimensionless radius and temperature graph.
function a9
R1=1;
R2=5;
eta=linspace(0,1,11);
zeta=linspace(0,1,11);
theta=linspace(0,pi,11);
R=linspace(1,5,11);
%Defining x-y boundary
for i=1:11
 x(1,i)= R1*cos(theta(i));
 y(1,i)= R1*sin(theta(i));
 
 x(11,i)= R2*cos(theta(i));
 y(11,i)= R2*sin(theta(i));
 
end
%Defining temperature boundary
t=zeros(11,11);
t(:,1)=30;
t(:,11)=100;
%Defining x-y boundary
for j=1:11
 x(j,1)= R(j)*cos(0);
 y(j,1)= 0;
 x(j,11)= R(j)*cos(pi);
 y(j,11)= 0;
end
%Interpolation Inside x and y

for i=2:10
    scalar_x=x(11,i)-x(1,i);
    scalar_y=y(11,i)-y(1,i);
    for j=2:10
      x(j,i)=x(1,i)+zeta(j)*scalar_x;  
      y(j,i)=y(1,i)+zeta(j)*scalar_y;  
    end
end
E1=1;
E2=1;
E3=1;


%To get a,b,c,x,y
delta_eta=eta(2)-eta(1);
delta_zeta=zeta(2)-zeta(1);
k=0;

while E1>10^-6 || E2>10^-6 || E3>10^-6
    
for j=2:10
    for i=2:10
x_eta(j-1,i-1)=(x(j,i+1)-x(j,i-1))/(2*delta_eta);
x_zeta(j-1,i-1)=(x(j+1,i)-x(j-1,i))/(2*delta_zeta);
y_eta(j-1,i-1)=(y(j,i+1)-y(j,i-1))/(2*delta_eta);
y_zeta(j-1,i-1)=(y(j+1,i)-y(j-1,i))/(2*delta_zeta);
J(j-1,i-1)=1/(x_eta(j-1,i-1)*y_zeta(j-1,i-1)-x_zeta(j-1,i-1)*y_eta(j-1,i-1));
a(j-1,i-1)=(J(j-1,i-1)*y_zeta(j-1,i-1)).^2+(-J(j-1,i-1)*x_zeta(j-1,i-1)).^2;
b(j-1,i-1)=-(J(j-1,i-1).^2)*(y_eta(j-1,i-1)*y_zeta(j-1,i-1)-x_eta(j-1,i-1)*x_zeta(j-1,i-1));
c(j-1,i-1)=(x_eta(j-1,i-1).^2+y_eta(j-1,i-1).^2)*J(j-1,i-1)^2;

X(j,i)=((a(j-1,i-1)/(delta_eta.^2))*(x(j,i+1)+x(j,i-1))+(b(j-1,i-1)/(2*delta_eta*delta_zeta))*(x(j+1,i+1)-x(j-1,i+1)-x(j+1,i-1)+x(j-1,i-1))+(c(j-1,i-1)/delta_zeta.^2)*(x(j+1,i)+x(j-1,i)))/(2*a(j-1,i-1)/(delta_eta.^2)+2*c(j-1,i-1)/(delta_zeta.^2));
Y(j,i)=((a(j-1,i-1)/(delta_eta.^2))*(y(j,i+1)+y(j,i-1))+(b(j-1,i-1)/(2*delta_eta*delta_zeta))*(y(j+1,i+1)-y(j-1,i+1)-y(j+1,i-1)+y(j-1,i-1))+(c(j-1,i-1)/delta_zeta.^2)*(y(j+1,i)+y(j-1,i)))/(2*a(j-1,i-1)/(delta_eta.^2)+2*c(j-1,i-1)/(delta_zeta.^2));


    end
end
%Temperature flux boundary condition and temperature calculations
T=zeros(11,11);
T(:,1)=t(:,1);
T(:,end)=t(:,end);
for i=2:10
  T(2,i)=((a(1,i-1)/(delta_eta.^2))*(t(2,i+1)+t(2,i-1))+(b(1,i-1)/(2*delta_eta*delta_zeta))*(t(3,i+1)-t(1,i+1)-t(3,i-1)+t(1,i-1))+(c(1,i-1)/delta_zeta.^2)*(t(3,i)+t(1,i)))/(2*a(1,i-1)/(delta_eta.^2)+2*c(1,i-1)/(delta_zeta.^2));  
end
T(1,:)=T(2,:);
for j=3:11
   T(j,:)=T(2,:); 
end
%Error calculation
err_t=T-t;
E3=max(max(abs(err_t)));
t=T;
xx=x(2:end-1,2:end-1);
XX=X(2:end,2:end);
x_r=x(1:end-1,1:end-1);
x_rr=x(1:end,1:end-1);
x_rrr=x(1:end,1:end-1);

 err_x=XX-xx;
 E1=max(max(abs(err_x)));
X_r=X;
X_r(1,:)=x_r(1,:);
X_r(:,1)=x_r(:,1);
X_r(11,:)=x_rr(end,:);
X_r(:,11)=x(:,end);

x=X_r;


yy=y(2:end-1,2:end-1);
YY=Y(2:end,2:end);
y_r=y(1:end-1,1:end-1);
y_rr=y(1:end,1:end-1);
y_rrr=y(1:end,1:end-1);

 err_y=YY-yy;
 E2=max(max(abs(err_y)));
Y_r=Y;
Y_r(1,:)=y_r(1,:);
Y_r(:,1)=y_r(:,1);
Y_r(11,:)=y_rr(end,:);
Y_r(:,11)=y(:,end);

y=Y_r;
k=k+1;
end
%Plotting
r=sqrt(x.^2+y.^2);
r_star=(r-R1)./(R2-R1);
t_star=(t'-30)./70;
%plot(r_star,t_star,'x');
plot(x,y,'r',x',y','r',R1*cos(theta),R1*sin(theta),'r',R2*cos(theta),R2*sin(theta),'r');
title('Grid generation');
figure
%Analytical solution plotting
T_1star=0.621334934*log(r);
plot(r_star,T_1star,'r');
title('Analytical Solution');
xlabel('rstar');
ylabel('Tstar');
figure
surf(x,y,T_1star);
 title('Analytical Solution-Body fitted Grid');
 xlabel('x');
 ylabel('y');
view (2)
colorbar
%plot(r_star,t_star-T_1star,'r');
 figure
 surf(x,y,t')
 view (2)
 title('Numerical Solution-Body fitted Grid');
 xlabel('x');
 ylabel('y');
 colorbar
end