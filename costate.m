function dlam=costate(t,y,S,E,I,U,t4)
beta=0.02763/1000;
rho=2.1;
tau=1/(60);
gamma=1/120;
k=0.2;
mu=1/(60.5);
n=0.9;
A=5000000;
S=interp1(t4,S,t);
E=interp1(t4,E,t);
I=interp1(t4,I,t);
U=interp1(t4,U,t);
dlam=zeros(4,1);
dlam(1)=y(1)*((1-U)*beta*(E+I)+ mu) +y(2)*(1-U)*beta*(E+I);
dlam(2)=-y(1)*(1-U)*beta*S-y(2)*((1-U)*beta*S-k*rho-tau*(1-k)-mu)+y(3)*k*rho...
    +y(4)*(tau*(1-k));
dlam(3)=-A-y(1)*(1-U)*beta*S-y(2)*(1-U)*beta*S+y(3)*(gamma+mu)+y(4)*gamma;
dlam(4)=y(4)*n*mu;

