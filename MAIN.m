%% Main function: integrating the dynamics equation of an elastic microfilament

%  Returns N+2 parameters for each time : 
%   x and y are the coordinates of the end of the first link
%   theta is the orientation of the first link
%   alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th
%   link

clear all;

%% number of links
N=40;

%% adimensionalizing parameters
global Hp Ho gamma Sp

Sp_c=4;
Sp=(Sp_c/N)^(4/3); %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta

% magnetic fields
Hx=@(t) 0;%100;
Hy=@(t) 0;%1000*cos(t/10);
Hp=@(z,t) 0;%cos(z(3))*Hx(t)+sin(z(3))*Hy(t);
Ho=@(z,t) 0;%-sin(z(3))*Hx(t)+cos(z(3))*Hy(t);

%% initial condition

% straight line
z0=zeros(N+2,1);

% half-circle
% z0=[N/(2*pi);0;-pi/2-pi/(2*N);-pi/(N)*ones(N-1,1)];

% quarter of a circle
% TT=[pi/(2*(N+1)):pi/(2*(N+1)):pi/2-pi/(2*(N+1))];
% x0=N*cos(TT-pi/(2*(N+1)));
% y0=N*sin(TT-pi/(2*(N+1)));
% t0=pi/2+TT;
% z0=[x0(1),y0(1),t0(1),t0(2:end)-t0(1:end-1)];

% parabola arc
% [x,y,th]=coord_parabola(N,-2,2);
% z0=[x(1),y(1),th(1),th(2:end)-th(1:end-1)];

%% time step
% T=1e-2;
% Nt=5000;
% dt=T/Nt;
T=8*pi; %final time
tps=linspace(0,T,800); %time step

%% Solving procedure
dZ=@(t,z) second_member_Nparam(t,z,N);
%options=odeset('RelTol',1e-2);
tic
[tps,traj]=ode15s(dZ,tps,z0);
toc

%% viz

% swimmer plot
figure(1)
hold on
col=colormap(jet(106));
A=[];
for i=length(tps)-150:15:length(tps)-45
    [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
    plot(X(1:end-1),Y(1:end-1),'--','LineWidth',2,'Color',col(i-(length(tps)-151),:))
    hold on
    % title(['N=',num2str(N)]);
    % axis equal
    axis([-0.1 1 -0.3 0.3])
    axis off
    drawnow
    % M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
end

figure(2)
clf
for i=400:length(tps)
    [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
    A(:,i)=TH(1:end);
end
surf(tps(end-400:end)/pi,[1:N]/N,A(:,end-400:end))
colormap jet
axis tight
shading interp
colorbar
% title('centerline angle')
xlabel('t / \pi')
ylabel('s')
set(gca,'FontSize',28)
% subplot(2,3,3)
% surf(tps(800:end),1:N-1,A)
% axis tight
% shading interp
% box on
% xlabel('time')
% ylabel('link number')
% colorbar 
% title('Curvature along the filament')
% subplot(2,3,2)
% plot(tps,TC)
% title('Total curvature')
% xlabel('time')
% ylabel('TC')
% 
% subplot(2,3,4)
% for i=1:30:300
%     [X,Y,TH]=coordinates_swimmerN(traj(i,:),N);
%     plot(X,Y,'Color',i/1.2/300*[1 1 1],'LineWidth',1)
%     hold on
% %     plot(X,Y,'Color',[0 0 0],'LineWidth',3)
% %     hold on
% %     plot(X1,Y1)
% %     plot(Xe,Ye,'LineWidth',3)
%     axis([-0.5 0.5 -0.1 0.5])
%     drawnow
% end
% subplot(2,3,6)
% surf(tps(1:300),1:N-1,A(:,1:300))
% axis tight
% shading interp
% box on
% xlabel('time')
% ylabel('link number')
% colorbar
% subplot(2,3,5)
% plot(tps(1:300),TC(1:300))
% xlabel('time')
% ylabel('TC')
