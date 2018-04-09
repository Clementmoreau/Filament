function B=second_member_Nparam_counterbend(t,z,N)
clear sum_th
th=zeros(N,1);
th(2)=z(3);

 for i=3:N+1
     th(i)=th(i-1)+z(i+1);
 end

B1=zeros(N+2,1);
amp=-0.4362/2;
a0p=amp*cos(t);
a0=amp*sin(t);
%bundle constant
KK=0.06;
%B1(N+2)=(th(N+1)-th(N))-(th(N)-a0)*KK;
B1(N+2)=(th(N+1)-th(N))+(th(N+1)-a0)*KK;
sum_th=th(N+1);
for i=N-2:-1:0
    sum_th=sum_th+th(i+2);
    %B1(3+i)=(th(i+2)-th(i+1))-(th(i+1)-a0)*KK;
    B1(3+i)=(th(i+2)-th(i+1))+(sum_th-(N-i)*a0)*KK;
end
 B1(3)=a0p;
% B1(3)=0;

BB=B1;
M=matrixNparam_counterbend(t,z,N);
B=M\BB;

end

