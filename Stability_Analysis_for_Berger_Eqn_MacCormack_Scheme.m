%This code plots range of values of D & C for stability for MacCormack
%Scheme. X-axis shows C values while Y- axis shows D values. MacCormack
%scheme is used to solve linear Burger's equation. 
%Berger's equation is unsteady partial differential equation  This numerical scheme is
%not stable for all the input values. This code helps to pick input values of C &
%D so that the scheme is stable.
function a6q3
theta=linspace(0,2*pi,20);
D=0.1;
C=0.1;
j=0;
 for D=0.1:0.01:1
      d=round(D*100);
 for C=0.1:0.01:1
      c=round(C*100);
   
for i=1:20
    A(i)=D^2*cos(2*theta(i))+(2*D+C^2-4*D^2)*cos(theta(i))+1-2*D-C^2+3*D^2;
    B(i)=(-C*D)*sin(2*theta(i))+(2*D-C)*sin(theta(i));
    g(i)=sqrt(A(i)^2+B(i)^2);
end
G(d,c)=max(abs(g));

if G(d,c)<1.000001
   j=j+1;
   dd(j)=d; 
   cc(j)=c;
   gg(j)=G(d,c);
   
end
 end
 end

GG=zeros(100,100);
 for h=1:length(dd)
GG(dd(h),cc(h))=1;
 end

  c=linspace(0.1,1,100);
  d=linspace(0.1,1,100);
  contour(c,d,GG)

 end