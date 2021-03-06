function cn(t1,t2,h,y0) 
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
k=ones(n+1,1);

for i=1:n
    x(i+1,1)=x(i,1)+h;
    ys=y(i,:);	
    B=[1-0.5*h*(-d-beta*ys(1,3)+a*(1-2*(ys(1,1)/Tmax))),-0.5*h*ro,0.5*h*beta*ys(1,1);-0.5*h*beta*ys(1,3),1+0.5*h*ro+0.5*h*delta,-0.5*h*beta*ys(1,1);0,-0.5*h*q,1+0.5*h*c];
    A=[ys(1,1)-y(i,1)-(h/2)*(f{1}(x(i,1),y(i,1),y(i,2),y(i,3))+f{1}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)));ys(1,2)-y(i,2)-(h/2)*(f{2}(x(i,1),y(i,1),y(i,2),y(i,3))+f{2}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)));ys(1,3)-y(i,3)-(h/2)*(f{3}(x(i,1),y(i,1),y(i,2),y(i,3))+f{3}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)))];
    yn=ys(1,:)-(B\A)';
    while ((norm(yn-ys)/norm(yn) > 10^-3) && k(i+1)<10)
        ys=yn;
        B=[1-0.5*h*(-d-beta*ys(1,3)+a*(1-2*(ys(1,1)/Tmax))),-0.5*h*ro,0.5*h*beta*ys(1,1);-0.5*h*beta*ys(1,3),1+0.5*h*ro+0.5*h*delta,-0.5*h*beta*ys(1,1);0,-0.5*h*q,1+0.5*h*c];
        A=[ys(1,1)-y(i,1)-(h/2)*(f{1}(x(i,1),y(i,1),y(i,2),y(i,3))+f{1}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)));ys(1,2)-y(i,2)-(h/2)*(f{2}(x(i,1),y(i,1),y(i,2),y(i,3))+f{2}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)));ys(1,3)-y(i,3)-(h/2)*(f{3}(x(i,1),y(i,1),y(i,2),y(i,3))+f{3}(x(i+1,1),ys(1,1),ys(1,2),ys(1,3)))];
        yn=ys(1,:)-(B\A)';
        k(i+1)=k(i+1)+1;
    end
    y(i+1,:)=yn;    
end
max(k);
end