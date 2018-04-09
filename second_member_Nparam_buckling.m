function B=second_member_Nparam_buckling(t,z,N,Sp)

th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+6,1);

B1(N+2)=(th(N+1)-th(N));

for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
end
B1(3)=0;
B1(N+3)=1.5; %buckling speed at proximal end
B1(N+5)=-1.5; %buckling speed at distal end


BB=B1;  
M=matrixNparam_buckling2(t,z,N,Sp);
B=M\BB;

end

