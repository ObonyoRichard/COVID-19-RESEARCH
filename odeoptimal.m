%Main function
%Authour obonyo Richard
%reference 
clear
clc
T=400;
tf=linspace(0,T,T+1);
U=zeros(1,T+1);
tspan1=tf;
IC=[1500 0 0.033 0];
a=500000;
beta1=0.02763/1000;
err=1;
umax=0.8;
%we solve our system with out control applied using 
%4th order runge-kutta
[T Y] = ode45(@odewithoutcontrol,tspan1,[1500 0 0.033 0]);
S1=Y(:,1);
E1=Y(:,2);
I1=Y(:,3);
R1=Y(:,4);

while err>0.00001;
    %we solve our state equations with control applied using 
%forward 4th order runge-kutta ode45
[T Y] = ode45(@(t,y)odewithcontrol(t,y,tf,U),tspan1,[1500 0 0.033 0]);
S=Y(:,1);
E=Y(:,2);
I=Y(:,3);
R=Y(:,4);
oldU=U;oldE=E;oldI=I;
%we flip the parameters up to down so as to use them to perform backward
%4th order runge-kutta using ode45.
c=flipud(S);
d=flipud(E);
e=flipud(I);
f=flipud(R);
v=flipud(U);
tspan2=flipud(T);
  %we solve our co-state equations with control applied using 
%backward 4th order runge-kutta. flipping the data so we have it all
%arranged from final time(tf=1001) to start time (t0=0)
[T Y] =  ode45(@(t,y)costate(t,y,c,d,e,v,tf),tspan2,[0 0 0 0]);
%we flip the values as the arary contain values from tfinal to initial
%time.
%inorder to compute our control to be applied we to re arrange the data 
%from starting time to end time
l1=flipud(Y(:,1));
l2=flipud(Y(:,2));
l3=flipud(Y(:,3));
l4=flipud(Y(:,4));
%we compute values of our control parameters to be applied for each time
%period.
for i=1:length(U)
    h=min((l1(i)-l2(i))*beta1*(E(i)+I(i))*S(i)/a,umax);
    Uval=max(0,h);
    U(i)=0.5*(Uval+oldU(i));
end
err=max(sum(abs(U))-sum(abs(oldU)),sum(abs(E))-sum(abs(oldE)));
end

infections_avoided=(sum(abs(E))-sum(abs(E1)))/sum(abs(E))

plot(tf,E1,'r',tf,E,'b')
legend('Without control','with control')
title('plot to show the effectiveness of optimal control')
ylabel('Number of infected(1000)')
xlabel('Days')
figure
plot(tf,I1,'r',tf,I,'b')
legend('Without control','with control')
ylabel('Number of infected(1000)')
xlabel('Days')
figure
plot(tf,U)
title('Agraph of optimal control parameter with time')
ylabel('Value of control applied')
xlabel('Days')