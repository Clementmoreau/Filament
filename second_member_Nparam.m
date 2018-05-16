% This function computes the value of \dot{X} at time t by computing A, Q and B from Eq. (17)

function B=second_member_Nparam(t,z,N)

th=zeros(N,1);
th(2)=z(3);

 for i=3:N+1
     th(i)=th(i-1)+z(i+1);
 end

B1=zeros(N+2,1);

B1(N+2)=(th(N+1)-th(N));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
end

 B1(3)=a0p;
% B1(3)=0;

BB=B1;
M=matrixNparam(t,z,N);
B=M\BB;

end

