% Main function: integrating the dynamics equation of a N-link elastic microswimmer
% Returns N+2 parameters : 
%   x and y are the coordinates of the end of the first link
%   theta is the orientation of the first link
%   alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th
%   link
function [tps,traj]=solver_swimmer_adaptative(N,Spe,T,delta_t,bar)

% adimensionalizing parameters
global gamma Sp 
Sp=Spe; %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta

%initial condition
z0f=[N/(2*pi);0;pi/2+pi/(2*N);pi/(N)*ones(N-1,1)];
[X0,Y0,TH0]=coordinates_swimmerN(z0f,N);
b=[sum(X0),sum(Y0)];
z0=[N/(2*pi)-(b(1)-bar(1));0-(b(2)-bar(2));pi/2+pi/(2*N);pi/(N)*ones(N-1,1)];
% z0=[N/pi;0;(pi/2+pi/(2*N));pi/N*ones(N-1,1)];
% z0=zeros(N+2,1);
% [x,y,th]=coord_parabola(N,-2,2);
% z0=[x(1),y(1),th(1),th(2:end)-th(1:end-1)];

%time step
% nt=T/delta_t;
% tps=[0:nt]*delta_t;
tps=linspace(0,T,1001);

dZ=@(t,z) second_member_Nparam(t,z,N);

%Solving procedure solver used could be ode1 ; ode2; ode3; ode5; ode15s
%options=odeset('RelTol',1e-1,'AbsTol',1e-4);
%tic
disp(['N = ',num2str(N)]);
[tps,traj]=ode15s(dZ,tps,z0);
%toc
end

