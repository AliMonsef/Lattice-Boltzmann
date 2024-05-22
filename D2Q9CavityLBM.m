%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Cavity Driven Flow With Lattice Boltzmann Method        %
%        This Program has been developed by "Ali Monsef"         %
%               Master of Mechanical Engineering                 %
%               Contact: alimonsef1997@gmail.com                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
%------------simulation parameters-------------
L=1;  U=1; dx=0.025;  nu=0.01;	% Physical kinematic viscosity
Re=L*U/nu;                      % physical Re number
NX=1/dx;   NY=NX;               % cavity length and cavity width
tau=0.7;   omega=1/tau;	        % relaxation time (BGK model) and omega

dt=((tau-0.5)*dx*dx)/(3*nu);	% time step
Ustar=U*dt/dx;               	% maximum velocity
nuL=(tau-0.5)/3;                % kinematic viscosity in Lattice space
ReL=NX*Ustar/nuL;               % Reynolds number in Lattice space
%------------Lattice parameters-------------
w  = [1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9]; % weights
cx = [1 0 -1  0 1 -1 -1  1 0];                  % velocities, x components
cy = [0 1  0 -1 1  1 -1 -1 0];                  % velocities, y components
% Node locations
x = (1:NX)-0.5; %link wise
y = (1:NY)-0.5; %link wise
% initialize populations
feq=zeros(NX,NY,9);
u=zeros(NX,NY);
v=zeros(NX,NY);
u_old=zeros(NX,NY);
rho=ones(NX,NY);
for k=1:9
    feq(:,:,k)=w(k);            % assuming density equal one and zero velocity initial state
end
f=feq;
fprop=feq;
tol=1e-5;                       % tolerance to steady state convergence
%-------------Main algorithm-------------
for t=1:15000
    for j=1:NY
        for i=1:NX
            for k=1:9
                feq(i,j,k)=rho(i,j)*w(k)*(1.0+3*(u(i,j)*cx(k)+v(i,j)*cy(k))+4.5*((u(i,j)*cx(k)+v(i,j)*cy(k))*(u(i,j)*cx(k)+v(i,j)*cy(k)))-1.5*(u(i,j)*u(i,j)+v(i,j)*v(i,j)));
            end
        end
    end
    f = (1-omega)*fprop + omega*feq;
    for k=1:9
        for j=1:NY
            for i=1:NX
                % Streaming step
                newx=1+mod(i-1+cx(k)+NX,NX);
                newy=1+mod(j-1+cy(k)+NY,NY);
                fprop(newx,newy,k)=f(i,j,k);
            end
        end
    end
    % Boundary condition (bounce-back)
    % Top wall (moving with tangential velocity Ustar)
    fprop(:,NY,4)=f(:,NY,2);
    fprop(:,NY,7)=f(:,NY,5)-(Ustar/6);
    fprop(:,NY,8)=f(:,NY,6)+(Ustar/6);
    % Bottom wall (rest)
    fprop(:,1,2)=f(:,1,4);
    fprop(:,1,5)=f(:,1,7);
    fprop(:,1,6)=f(:,1,8);
    % Right wall (rest)
    fprop(NX,:,3)=f(NX,:,1);
    fprop(NX,:,6)=f(NX,:,8);
    fprop(NX,:,7)=f(NX,:,5);
    % Left wall (rest)
    fprop(1,:,1)=f(1,:,3);
    fprop(1,:,5)=f(1,:,7);
    fprop(1,:,8)=f(1,:,6);
    % Compute macroscopic quantities
    for i=1:NX
        for j=1:NY
            rho(i,j) =fprop(i,j,1)+fprop(i,j,2)+fprop(i,j,3)+fprop(i,j,4)+fprop(i,j,5)+fprop(i,j,6)+fprop(i,j,7)+fprop(i,j,8)+fprop(i,j,9);
            u(i,j) =(fprop(i,j,1)+fprop(i,j,5)+fprop(i,j,8)-(fprop(i,j,3)+fprop(i,j,6)+fprop(i,j,7)))/rho(i,j);
            v(i,j) =(fprop(i,j,2)+fprop(i,j,5)+fprop(i,j,6)-(fprop(i,j,4)+fprop(i,j,7)+fprop(i,j,8)))/rho(i,j);
        end
    end
    % check convergence
    if mod(t,200)==1
        t;
        Err = abs(u-u_old);
        if norm(Err,1)<tol
            break
        else
            u_old = u;
        end
    end
end
%------------plot-------------
u=rot90(u,3);
u=fliplr(u);
v=rot90(v,3);
v=fliplr(v);
figure(1);streamslice(x,y,u,v)
title('Stream Lines')
figure(2); contour(x,y,u)
title('u velocity contour')
