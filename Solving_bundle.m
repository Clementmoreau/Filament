% Main function: integrating the dynamics equation of a N-link elastic bundle made of two filaments

% number of links
N=40;

% adimensionalizing parameters
global gamma Sp k a th_max
a=1;
th_max=0.28*pi;

Sp_c=4;
Sp=(Sp_c/N)^(4/3); %sperm number characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/3) 
gamma=1/2; %ratio between the hydrodynamics drag coefficients - gamma = xi/eta
k=1*N/(a^2)/1000; %ratio between the spring stiffnesses

% initial condition
z0=[[0:N-1]';zeros(N,1);zeros(N,1);[0:N-1]';a*ones(N,1);zeros(N,1)];
% TT=[pi/(2*(N+1)):pi/(2*(N+1)):pi/2-pi/(2*(N+1))];
% x0=N*cos(TT-pi/(2*(N+1)));
% y0=N*sin(TT-pi/(2*(N+1)));
% t0=pi/2+TT;
% z0=[x0';y0';t0';x0';y0'+a;t0'];

% time step
T=8*pi;%resolution time
tps=linspace(0,T,800);

% Solving procedure
dZ=@(t,z) second_member_bundle(t,z,N);
tic
[tps,traj]=ode15s(dZ,tps,z0);
toc

% viz
figure(1)
clf
hold on
col=colormap(jet(106));
for i=length(tps)-150:15:length(tps)-45

plot(([traj(i,1:N)]/N+[traj(i,3*N+1:4*N)]/N)/2,([traj(i,N+1:2*N)]/N+[traj(i,4*N+1:5*N)]/N)/2-a/2/N,...
        'Color',col(i-(length(tps)-151),:),'LineWidth',2)
end

%%

% swimmer plot
figure(2)
clf
col=colormap(jet(106));
for i=length(tps)-150:15:length(tps)-45
    plot([traj(i,1:N)]/N,[traj(i,N+1:2*N)]/N,...
        'Color',col(i-(length(tps)-151),:),'LineWidth',3)
    hold on
    plot([traj(i,3*N+1:4*N)]/N,[traj(i,4*N+1:5*N)]/N,...
        'Color',col(i-(length(tps)-151),:),'LineWidth',3)
    for j=1:N-1
        plot([traj(i,j+1),traj(i,3*N+j+1)]/N,[traj(i,N+j+1),traj(i,4*N+j+1)]/N,...
        'Color',col(i-(length(tps)-151),:))
    end
    axis([-0.1 1 -0.3 0.3])
    axis off
    drawnow
    % M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.

end

nb_per=200*2;

% centerline angle

figure(3)
clf
surf(tps(end-nb_per:end)/pi,[1:N]/N,(traj(end-nb_per:end,2*N+1:3*N)'+traj(end-nb_per:end,5*N+1:end)')/2)
colormap jet
axis tight
shading interp
colorbar
% title('centerline angle')
xlabel('t / \pi')
ylabel('s')
set(gca,'FontSize',28)

%sliding displacement 


sldx = (traj(end-nb_per:end,1:N)-traj(end-nb_per:end,3*N+1:4*N))/N;
sldy = (traj(end-nb_per:end,N+1:2*N)-traj(end-nb_per:end,4*N+1:5*N))/N;
sld = sqrt(sldx.^2+sldy.^2);
figure(4)
clf
surf(tps(end-nb_per:end)/pi,[1:N]/N,(sld(:,1:N)'-a/N)*N)
colormap(jet(2000))
axis tight
shading interp
colorbar
caxis([-0.3 0.2])
% title('sliding displacement (difference btw length btw 2 nodes and the rest length, as a percentage of the rest length)')
xlabel('t / \pi')
ylabel('s')
set(gca,'FontSize',28)


%distance between filaments
d1x = (traj(end-nb_per:end,2:N)-traj(end-nb_per:end,1:N-1))/N;
d1y = (traj(end-nb_per:end,N+2:2*N)-traj(end-nb_per:end,N+1:2*N-1))/N;
dbf = N*abs(sldx(:,1:end-1).*d1y-sldy(:,1:end-1).*d1x);
figure(5)
clf
surf(tps(end-nb_per:end)/pi,[1:N-1]/(N-1),N*dbf'-1)
colormap(jet(2000))
axis tight
shading interp
colorbar
caxis([-0.3 0.2])
% title('dist btw filaments as a percentage of the rest distance')
xlabel('t / \pi')
ylabel('s')
set(gca,'FontSize',28)





