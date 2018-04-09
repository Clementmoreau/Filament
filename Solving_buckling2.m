%% number of links
N=20;

%% adimensionalizing parameters
global gamma Sp th0 f0

Sp_c=4;
Sp=(Sp_c/N)^(4/3); %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta
f0=1e-2;
th0=pi/6;

%% initial condition
z0=[0;-1;th0;-2e-2*ones(N-1,1)];

%% time step
% T=1e-2;
% Nt=5000;
% dt=T/Nt;
T=3;%resolution time
tps=linspace(0,T,300);

%% Solving procedure
dZ=@(t,z) second_member_Nparam_buckling2(t,z,N);
%options=odeset('RelTol',1e-2);
tic
[tps,traj]=ode15s(dZ,tps,z0);
toc

%% viz
figure;
for i=1:length(tps)
    [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
    plot(X,Y,'k','LineWidth',2)
    axis equal
    title(['T = ',num2str(tps(i))]);
    drawnow
end