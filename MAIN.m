% *******************************************
% Main function: integrating the dynamics equation of an elastic microfilament
% *******************************************

%  References to the original paper are given throughout the code
%
%  Returns N+2 parameters for each time (contained in the 'traj' array) : 
%   x and y are the coordinates of the end of the first link
%   theta is the orientation of the first link
%   alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th link
%  See Fig.1 

% -----------------
clear all;

% ---- Choose number of links
N=40;

% Adimensionalizing parameters
global Hp Ho gamma Sp

% ---- Choose 'sperm number'characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/4)
Sp_c=4;

Sp=(Sp_c/N)^(4/3);
gamma=1/2; % ratio between the hydrodynamics drag coefficients - gamma = xi/eta

% ---- Choose magnetic fields if relevant (see section V-B; example in comments)
Hx=@(t) 0; %100;
Hy=@(t) 0; %1000*cos(t/10);
Hp=@(z,t) cos(z(3))*Hx(t)+sin(z(3))*Hy(t);
Ho=@(z,t) -sin(z(3))*Hx(t)+cos(z(3))*Hy(t);

% ---- Choose initial condition (uncomment the require one)

% Straight line
z0=zeros(N+2,1);

% Half-circle
% z0=[N/(2*pi);0;-pi/2-pi/(2*N);-pi/(N)*ones(N-1,1)];

% Parabola arc
% [x,y,th]=coord_parabola(N,-2,2);
% z0=[x(1),y(1),th(1),th(2:end)-th(1:end-1)];

% ---- Choose final time and time step
T=8*pi; %final time
tps=linspace(0,T,800); %time step

% ---- Solving procedure
% ---- Choose a function 
%   * second_member_Nparam is for the standard nonmagnetised case (section IV)
%   *
%   *
dZ=@(t,z) second_member_Nparam(t,z,N);
% tic
[tps,traj]=ode15s(dZ,tps,z0);
% toc

% ---- Graphic visualisation 
figure;


