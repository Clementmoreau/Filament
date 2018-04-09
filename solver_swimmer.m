%% Main function: integrating the dynamics equation of a N-link elastic microswimmer
% Returns N+2 parameters : 
%   x and y are the coordinates of the end of the first link
%   theta is the orientation of the first link
%   alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th
%   link
function [tps,traj]=solver_swimmer(N,Spe,T,delta_t)

% adimensionalizing parameters
global Hx Hy gamma Sp %LL
Sp=Spe; %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta
%LL=2*sin(pi/2/N);
%magnetic fields
Hx=@(t) 0;%100;
Hy=@(t) 0;%1000*cos(t/10);

%initial condition
z0=[N/pi;0;pi/2+pi/(2*N);pi/N*ones(N-1,1)];
%z0=zeros(N+2,1);

%time step
% T=1e-2;
% Nt=5000;
% dt=T/Nt;

dZ=@(t,z) second_member_Nparam(t,z,N);
tps=linspace(0,T,T/delta_t);

%Solving procedure solver used could be ode1 ; ode2; ode3; ode5; ode15s
%options=odeset('RelTol',1e-1,'AbsTol',1e-4);
%tic
disp(['N = ',num2str(N)]);
traj=ode5(dZ,tps,z0);
%toc
