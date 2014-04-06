function [y]=RK4(t1,t2,h,y0)
%t1-pociatocny cas, t2-koncovy cas, h-krok, y0-zaciatocna podmienka(riadkovy vektor)
s=5;
d=0.01;
a=0.5;
Tmax=1200;
beta=0.0002;
ro=0.01;
delta=1;
q=800;
c=5;

f{1}=@(t,y1,y2,y3) s-d*y1+a*y1*(1-(y1/Tmax))-beta*y1*y3+ro*y2;
f{2}=@(t,y1,y2,y3) beta*y1*y3-y2*delta-y2*ro;
f{3}=@(t,y1,y2,y3) q*y2-c*y3;

n=(t2-t1)/h;
x=zeros(n+1,1);
x(1,1)=t1;
y=zeros(n+1,3);
y(1,:)=y0;

k1=zeros(1,3);
k2=zeros(1,3);
k3=zeros(1,3);
k4=zeros(1,3);

for i=1:n
    for j=1:3
        k1(j)=h*f{j}(x(i,1),y(i,1),y(i,2),y(i,3));
        k2(j)=h*f{j}(x(i,1)+0.5*h,y(i,1)+0.5*k1(j),y(i,2)+0.5*k1(j),y(i,3)+0.5*k1(j));
        k3(j)=h*f{j}(x(i,1)+0.5*h,y(i,1)+0.5*k2(j),y(i,2)+0.5*k2(j),y(i,3)+0.5*k2(j));
        k4(j)=h*f{j}(x(i,1)+h,y(i,1)+k3(j),y(i,2)+k3(j),y(i,3)+k3(j));
        y(i+1,j)=y(i,j)+(1/6)*(k1(j)+2*k2(j)+2*k3(j)+k4(j));
    end
    x(i+1,1)=x(i,1)+h;
end

end