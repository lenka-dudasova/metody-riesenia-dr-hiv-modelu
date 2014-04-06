function AdamsPK(t1,t2,h,y0)
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

kt=RK4(t1,t1+3*h,h,y0,p);

for i=2:4
	y(i,:)=kt(i,:);
    x(i+1,1)=x(i,1)+h;

end

for i=4:n
    yv=zeros(1,3);
    x(i+1,1)=x(i,1)+h;
    for j=1:3
        yv(1,j)= y(i,j)+(h/24)*(55*f{j}(x(i),y(i,1),y(i,2),y(i,3))-59*f{j}(x(i-1),y(i-1,1),y(i-1,2),y(i-1,3))+37*f{j}(x(i-2),y(i-2,1),y(i-2,2),y(i-2,3))-9*f{j}(x(i-3),y(i-3,1),y(i-3,2),y(i-3,3)));
    end	
	
    for k=1:2
		a1=f{1}(x(i),y(i,1),y(i,2),y(i,3));
		a2=f{1}(x(i-1),y(i-1,1),y(i-1,2),y(i-1,3));
		a3=f{1}(x(i-2),y(i-2,1),y(i-2,2),y(i-2,3));
		b1=f{2}(x(i),y(i,1),y(i,2),y(i,3));
		b2=f{2}(x(i-1),y(i-1,1),y(i-1,2),y(i-1,3));
		b3=f{2}(x(i-2),y(i-2,1),y(i-2,2),y(i-2,3)); 
		c1=f{3}(x(i),y(i,1),y(i,2),y(i,3));
		c2=f{3}(x(i-1),y(i-1,1),y(i-1,2),y(i-1,3));
		c3=f{3}(x(i-2),y(i-2,1),y(i-2,2),y(i-2,3));
		C=[((h/24)*(9*f{1}(x(i+1),yv(1,1),yv(1,2),yv(1,3))+19*a1-5*a2+a3)),((h/24)*(9*f{2}(x(i+1),yv(1,1),yv(1,2),yv(1,3))+19*b1-5*b2+b3)),((h/24)*(9*f{3}(x(i+1),yv(1,1),yv(1,2),yv(1,3))+19*c1-5*c2+c3))];
		yn(1,:)=y(i,:)+C(1,:);
		yv=yn;
    end
	y(i+1,:)=yn;    
end

end