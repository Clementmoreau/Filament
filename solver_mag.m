function [tps,traj]=solver_mag(N,Spe,T,delta_t,M)

global gamma Sp 

% initial straight line
z0=zeros(N+2,1);

%time step
nt=T/delta_t;
tps=[0:nt]*delta_t;

dZ=@(t,z) second_member_mag(t,z,N,M);

[tps,traj]=ode15s(dZ,tps,z0);

end