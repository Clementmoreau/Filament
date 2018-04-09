function B=second_member_Nparamfixe(t,z,N)
global L Mag Eta Xi k Hx Hy Sp
Mag=zeros(1,N+3);
th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B2=zeros(N+2,1);
B3=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));
B2(N+2)=-Mag(N)*sin(th(N));
B3(N+2)=Mag(N)*cos(th(N));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
    B2(3+i)=B2(3+i+1)-Mag(i+1)*sin(th(i+2));
    B3(3+i)=B3(3+i+1)+Mag(i+1)*cos(th(i+2));
end
B1(3)=2*sin(pi*t)*Sp;

BB=B1+Hx(t)*B2+Hy(t)*B3;
M=matrixNparamfixe(t,z,N);
B=M\BB;

end

