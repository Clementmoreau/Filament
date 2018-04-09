function B=second_member_Nparam_oscillation(t,z,N)
th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
end
B1(3)=-1*cos(3*(2*pi)*t/10)/N;
BB=B1;
M=matrixNparam_oscillation(t,z,N);
B=M\BB;

end

