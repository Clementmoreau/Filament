% This function solves the linear system (20)
% from Appendix VII-C in the case of a pinned end with a forced
% angular actuation

function B=second_member_Nparam_oscillation(t,z,N)
th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));

% for forced angular actuation
% --- Choose an amplitude amp in rad
amp=0.28*pi;
a0p=amp*cos(t);
a0=amp*sin(t);

for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
end

B1(3)=a0p;

M=matrixNparam_oscillation(t,z,N);
B=M\B1;

end

