% *******************************************
% Main function: integrating the dynamics equation of an elastic microfilament
% *******************************************
%  References to the original paper are given throughout the code
%
%  Returns N+2 parameters for each time (contained in the 'traj' array) : 
%   * x and y are the coordinates of the end of the first link
%   * theta is the orientation of the first link
%   * alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th link
%  See Fig.1 
% *******************************************
%
clear all;

% ---- Choose number of links
N=40;

% Declaring global parameters
global gamma Sp

% ---- Choose 'sperm number' characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/4)
% typically between 2 and 16, see Eq. (5)
Sp_c=4;
Sp=(Sp_c/N)^(4/3);
gamma=1/2; % ratio between the hydrodynamics drag coefficients - gamma = xi/eta

% ---- Choose initial condition (uncomment the required one)

% (1) Straight line
 z0=zeros(N+2,1);

% (2) Half-circle
% z0=[N/(2*pi);0;-pi/2-pi/(2*N);-pi/(N)*ones(N-1,1)];

% (3) Parabola arc
% [x,y,th]=coordinates_parabola(N,-2,2);
% z0=[x(1),y(1),th(1),th(2:end)-th(1:end-1)];

% ---- Choose final time and time step
T=8; %final time
tps=linspace(0,T,200); %time step

% ---- Choose a function 1, 2 or 3 (uncomment the required section)

% (1) relaxation is for the standard nonmagnetised case (section IV)
% in that case start from a non straight configuration
 dZ=@(t,z) relaxation(t,z,N);
% ***********

% (2) oscillation is for the standard case with a pinned end and angular actuation
% % ---- Choose an angular amplitude in rad
% amp=0.5*pi;
% dZ=@(t,z) oscillation(t,z,N,amp);
% ************

% (3) magnetism is for the magnetised case (section V-B)
% ---- Choose a magnetisation for the filament (N values)
% Mag=[ones(1,N/4),zeros(1,3*N/4)]; 
% % ---- Choose time-varying magnetic fields Hx and Hy (see section V-B)
% (I recommend not to use fields with absolute value greater than 1)
% Hx=@(t) 0;
% Hy=@(t) 0.03*cos(t*10);
% dZ=@(t,z) magnetism(t,z,N,Mag,Hx,Hy);
% ************

% ---- SOLVING
[tps,traj]=ode15s(dZ,tps,z0);

% ---- Graphic visualisation 
figure;
for i = 1:length(tps)
     [X,Y,TH]=coordinates_filament(traj(i,:),N);
     plot(X,Y,'k','LineWidth',0.5)
     axis([-1 1 -1 1])
     title(['T = ',num2str(tps(i))]);
     drawnow
end

% *********** END OF MAIN SCRIPT ***********
% ***** SUMMARY OF FOLLOWING FUNCTIONS *****
% relaxation
% oscillation
% magnetism
% matrixNparam
% matrixNparam_oscillation
% matrix3Nparameters
% coordinates_filament
% coordinates_parabola
% *******************************************

function B=relaxation(t,z,N)
% This function computes the value of \dot{X} at time t by computing A, Q and B from Eq. (20)
% see Appendix VII-C, Eq. (20)

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

B1(3)=0;

BB=B1;

% call to the function that computes Sp^4*A*Q
M=matrixNparam(t,z,N);

% solving the linear system
B=M\BB;

end

% **************************************
function B=oscillation(t,z,N,amp)
% This function solves the linear system (20)
% from Appendix VII-C in the case of a pinned end with a forced
% angular actuation

th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B1(N+2)=(th(N+1)-th(N));

% for forced angular actuation
a0p=amp*cos(t);

for i=N-2:-1:0
    B1(3+i)=(th(i+2)-th(i+1));
end

B1(3)=a0p;

M=matrixNparam_oscillation(t,z,N);
B=M\B1;

end

% **************************************
function B=magnetism(t,z,N,Mag,Hx,Hy)
% This function is similar to second_member, but with the magnetic effects added
% (see Appendix VII-E, Eq. 24)
% Hx and Hy are the external magnetic fields, defined in the MAIN file.

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

% **************************************
function M=matrixNparam(t,z,N)
% This function fills the matrix Q as defined in the text
% (see Appendix VII-B and C, equation (19) and (20))
% and returns the product Sp^4*A*Q

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+cos(z3(2*N+i-1));
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1));
end

M3=matrix3Nparameters(t,z3,N);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;
    end
end
for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1));
    C2(i,i-1)=cos(z3(2*N+i-1));
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j));
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j));
    end
end
C=[ones(N,1),zeros(N,1),C1;zeros(N,1),ones(N,1),C2;zeros(N,2),C3];
M=M3*C;

end

% **************************************
function  M=matrixNparam_oscillation(t,z,N)
% This function is similar to matrixNparam but in the 
% case of a pinned proximal end with a forced angular actuation.

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+cos(z3(2*N+i-1));
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1));
end

M3=matrix3Nparameters(t,z3,N);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;
    end
end
for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1));
    C2(i,i-1)=cos(z3(2*N+i-1));
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j));
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j));
    end
end
C=[ones(N,1),zeros(N,1),C1;zeros(N,1),ones(N,1),C2;zeros(N,2),C3];
M=M3(1:N+2,:)*C;

% --- this is where the first two equations are being changed to 
% implement the pinned proximal end

M(1,:)=[1,zeros(1,N+1)];
M(2,:)=[0,1,zeros(1,N)];
M(3,:)=[0,0,1,zeros(1,N-1)];
end

% **************************************
function M=matrix3Nparameters(t,z,N)
% This function fills the matrix defined as A in the text 
% (see Appendix VII-C, equation (20) and following)
% and returns Sp^4 * A

global gamma Sp

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);

for i=1:N
    u=cos(th(i));
    v=sin(th(i));
    F(1,i)=-(gamma*u^2+v^2);
    F(1,N+i)=-(gamma-1)*u*v;
    F(2,i)=-(gamma-1)*u*v;
    F(2,N+i)=-(u^2+gamma*v^2);
    F(1,2*N+i)=1/2*v;
    F(2,2*N+i)=-1/2*u;
end

F=Sp^3*F;

for i=1:N
    for j=i:N
        u=cos(th(j));
        v=sin(th(j));
        A=(x(j)-x(i));
        B=(y(j)-y(i));
        T(i,j)=1/2*v+...
            A*(-gamma+1)*v*u+...
            B*(gamma*u*u+v*v);
        T(i,N+j)=-1/2*u+...
            B*(gamma-1)*v*u-...
            A*(u*u+gamma*v*v);
        T(i,2*N+j)=-1/3-...
            1/2*A*u...
            -1/2*B*v;
    end
end

T=Sp^3*T;

M=[F;T];
end

% **************************************
function [X,Y,TH]=coordinates_filament(z,N)
% This function computes the '3N coordinates' -- X_3N in the text
% from the 'N+2 coordinates' -- X in the text 

% --- input : N+2 coordinates, number of links N
% --- output : X, Y coordinates of the N links
% TH orientation of each link

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

% **************************************
function [X,Y,TH]=coordinates_parabola(N,x1,x2)
% This function computes the coordinates of N links of constant length on a parabola y=x^2 arc
% between x=x1 and x=x2
% normalized at the end such that total length=N
% Used for the parabola initial condition

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
end
