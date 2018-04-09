%% Main function: integrating the dynamics equation of a N-link elastic microswimmer
% Returns N+2 parameters : 
%   x and y are the coordinates of the end of the first link
%   theta is the orientation of the first link
%   alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th
%   link

clear all;

%% number of links
N=36;

%% adimensionalizing parameters
global Hp Ho gamma Sp

Sp_c=2;
Sp=(Sp_c/N)^(4/3); %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta

%magnetic fields
Hx=@(t) 1e-3;%100;
Hy=@(t) 1e-2*cos(t);%1000*cos(t/10);
Hp=@(z,t) 0;
Ho=@(z,t) 0;

%% initial condition

% straight line
% z0=zeros(N+2,1);

% "fake" straight line for buckling
% p=randi([2 N-1]);
% s=rand;
% absbias=1e-3;
% bias=absbias*(s>1/2)-absbias*(s<=1/2);
% z0=[zeros(p,1);bias;-2*bias;bias;zeros(N-p-1,1);0;0;0;0];

absbias=1e-2;
% first mode (triangular)
 z0=[0;0;absbias;zeros(1,17)';-2*absbias;zeros(1,17)';zeros(4,1)];

% second mode
% z0=[0;0;absbias;zeros(1,8)';-2*absbias;zeros(1,17)';2*absbias;zeros(1,8)';zeros(4,1)];

% third mode
% z0=[0;0;absbias;zeros(1,5)';-2*absbias;zeros(1,11)';2*absbias;zeros(1,11)';-2*absbias;zeros(1,5)';zeros(4,1)];

%% time step
% T=1e-2;
% Nt=5000;
% dt=T/Nt;
T=1.4;%resolution time
tps=linspace(0,T,300);

%% Solving procedure
dZ=@(t,z) second_member_Nparam_buckling(t,z,N);
%options=odeset('RelTol',1e-2);
tic
[tps,traj]=ode15s(dZ,tps,z0);
toc

%% viz

% swimmer plot
Xe=[];
Ye=[];
A=[];
figure(1);
clf
col = colormap(jet(300));
hold off
%subplot(1,2,k)
for i=1:30:length(tps)
    [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
    %A(:,i)=TH(2:end)-TH(1:end-1);
    % plot(X,Y-(i-1)/800,'Color',i/1.4/length(tps)*[1 1 1],'LineWidth',3)
    plot(X,Y,'Color',col(i,:),'LineWidth',2)
    hold on
    axis([0 1 -0.2 0.2])
    axis off
end
figure(2)
clf
colormap(jet(1000))
for i=1:length(tps)
    A=[A;traj(i,4:N+2)];
end
s=surf(tps(1:end),[1:N-1]/(N-1),A');
caxis([-0.18 0.18]);
axis tight
colorbar off
c=colorbar;
c.Box='off';
c.TickLength=0;
shading interp
set(gca,'FontSize',30)
xlabel('t')
ylabel('s')

% figure(3)
% [X,Y,TH]=coordinates_swimmerN(traj(1,:),N);
% plot(X,Y,'Color','k','LineWidth',3)
% axis([0 1 -0.02 0.02])
% axis off



%% Movie
% figure;
% for i = 1:length(tps)
%     [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
%     % Example of plot
%     plot(X,Y,'k','LineWidth',0.5)
%     axis([0 1 -0.2 0.2])
%     title(['T = ',num2str(tps(i))]);
% 
%     % Store the frame
%     M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
% end
% 
% v=VideoWriter('BucklingSp2.avi');
% open(v);
% writeVideo(v,M);
% close(v);
