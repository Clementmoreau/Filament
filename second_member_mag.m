% This function is similar to second_member, but with the magnetic effects added
% (see Appendix VII-E, Eq. 24)
% Hx and Hy are the external magnetic fields, defined in the MAIN file.

function B=second_member_mag(t,z,N,Mag)
global Ho Hp Hx Hy

% --- Choose a magnetisation
Mag=[zeros(1,3*N/4),ones(1,N/4)];
th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B2=zeros(N+2,1);
B3=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));

% adding the magnetic effects 
B2(N+2)=-Mag(N)*sin(th(N));
B3(N+2)=Mag(N)*cos(th(N));
for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
    B2(3+i)=B2(3+i+1)-Mag(i+1)*sin(th(i+2));
    B3(3+i)=B3(3+i+1)+Mag(i+1)*cos(th(i+2));
end
B1(3)=0;

BB=B1+Hx(t)*B2+Hy(t)*B3;
M=matrixNparam(t,z,N);
B=M\BB;

end
