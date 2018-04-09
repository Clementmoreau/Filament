function B=second_member_Nparam_buckling(t,z,N,Sp)

th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+6,1);
% B2=zeros(N+2,1);
% B3=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));
% B2(N+2)=-Mag(N)*sin(th(N));
% B3(N+2)=Mag(N)*cos(th(N));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
%     B2(3+i)=B2(3+i+1)-Mag(i+1)*sin(th(i+2));
%     B3(3+i)=B3(3+i+1)+Mag(i+1)*cos(th(i+2));
end
B1(3)=0;
B1(N+3)=1.5;
B1(N+5)=-1.5;


BB=B1;   %+Hx(t)*B2+Hy(t)*B3;
M=matrixNparam_buckling2(t,z,N,Sp);
B=M\BB;

end

