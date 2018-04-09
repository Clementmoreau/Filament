function [X,Y,TH]=coordinates_swimmerN(z,N)

X=zeros(N+1,1);
Y=zeros(N+1,1);
TH=zeros(N,1);

X(1)=z(1);
Y(1)=z(2);
TH(1)=z(3);

for i=2:N
    X(i)=X(i-1)+cos(TH(i-1));
    Y(i)=Y(i-1)+sin(TH(i-1));
    TH(i)=TH(i-1)+z(i+2);
end

X(N+1)=X(N)+cos(TH(N));
Y(N+1)=Y(N)+sin(TH(N));
X=X/N;
Y=Y/N;

end
