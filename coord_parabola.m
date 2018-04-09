function [X,Y,TH]=coord_parabola(N,x1,x2)

%compute the coordinates of N links of constant length on a parabola y=x^2 arc
%between x=x1 and x=x2
%normalized at the end such that total length=N

F=@(x) (x*sqrt(1+x^2)+log(x+sqrt(1+x^2)))/2;
f=@(xa,xb) abs(xa-xb)*sqrt(1+(xa+xb)^2);

L_arc=(F(2*x2)-F(2*x1))/2;

% dichotomy initialization
a=L_arc/N;
b=L_arc/N/2;

x=x1;
tol=1e-10;

while abs(x-x2)>tol
    c=(a+b)/2;
    x=x1;
    %abscissa of the last point
    
    for i=1:N
        A=x;
        B=x+c;
        C=B;
        while abs(f(x,C)-c)>tol
            C=(A+B)/2;
            if (f(x,C)-c)>0
                B=C;
            else
                A=C;
            end
        end
        x=C;       
    end
    %dichotomy step
    if (x-x2) > 0
        a=c;
    else
        b=c;
    end
end

%coordinates
X(1)=x1;
x=x1;
for i=1:N
        A=x;
        B=x+c;
        C=B;
        while abs(f(x,C)-c)>tol
            C=(A+B)/2;
            if (f(x,C)-c)>0
                B=C;
            else
                A=C;
            end
        end
        X=[X C];
        x=C;
end

Y=X.^2;
TH=atan((Y(2:end)-Y(1:end-1))./(X(2:end)-X(1:end-1)));

l=f(X(1),X(2));
X=X/(l);
Y=Y/(l);


% figure;
% plot(X,Y)

end




    