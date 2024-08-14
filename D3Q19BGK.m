%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           D3Q19 (BGK) With Lattice Boltzmann Method            %
%        This Program has been developed by "Ali Monsef"         %
%               Master of Mechanical Engineering                 %
%               Contact: alimonsef1997@gmail.com                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
%------------simulation parameters-------------
nx = 15; % length
ny = nx; % width
nz = nx; % altitude
omega = 1.0;
density = 1.0;
%---weights for the D3Q19 lattice model---
t1 = 1/3; 
t2 = 1/18;
t3 = 1/36;
f = repmat(density/19,[nx ny nz 19]);
feq = f; % equilibrium distribution
matsize = nx * ny * nz;
CI = 0 : matsize : matsize * 19 ;
boundary = zeros(nx,ny,nz); % boundary conditions

%-------------Main algorithm-------------
for i=1:nx
    for j=1:ny
        for k=1:nz
        	boundary(i,j,k) = ((i-5)^2+(j-6)^2+(k-7)^2)<6;
        end
    end
end
boundary(:,:,1)=1;
boundary(:,1,:)=1;
on=find(boundary); %---matrix offset of each Occupied Node---

TO_REFLECT=[on+CI(2) on+CI(3) on+CI(4) on+CI(5)	on+CI(6) on+CI(7) on+CI(8) ...
    on+CI(9) on+CI(10) on+CI(11) on+CI(12) on+CI(13) on+CI(14) on+CI(15) ...
    on+CI(16) on+CI(17) on+CI(18) on+CI(19)];
reflected=[on+CI(3) on+CI(2) on+CI(5) on+CI(4) on+CI(7) on+CI(6) on+CI(11) ...
    on+CI(10) on+CI(9) on+CI(8) on+CI(15) on+CI(14) on+CI(13) on+CI(12) ...
    on+CI(19) on+CI(18) on+CI(17) on+CI(16)];
avu=1; prevavu=1; ts=0; deltaU=1e-7; numactivenodes=sum(sum(sum(1-boundary)));

while (ts<4000 & 1e-10<abs((prevavu-avu)/avu)) | ts<100
    %---Propagate---
    %---nearest-neighbours---
    f(:,:,:,2)=f(:,:,[nz 1:nz-1],2);
    f(:,:,:,3)=f(:,:,[2:nz 1],3);
    f(:,:,:,4)=f(:,[ny 1:ny-1],:,4);
    f(:,:,:,5)=f(:,[2:ny 1],:,5);
    f(:,:,:,6)=f([nx 1:nx-1],:,:,6);
    f(:,:,:,7)=f([2:nx 1],:,:,7);
    %---next-nearest neighbours---
    f(:,:,:,8)= f([nx 1:nx-1],[ny 1:ny-1],:,8);
    f(:,:,:,9)= f([nx 1:nx-1],[2:ny 1],:,9);
    f(:,:,:,10)=f([2:nx 1],[ny 1:ny-1],:,10);
    f(:,:,:,11)=f([2:nx 1],[2:ny 1],:,11);
    f(:,:,:,12)=f([nx 1:nx-1],:,[nz 1:nz-1],12);
    f(:,:,:,13)=f([nx 1:nx-1],:,[2:nz 1],13);
    f(:,:,:,14)=f([2:nx 1],:,[nz 1:nz-1],14);
    f(:,:,:,15)=f([2:nx 1],:,[2:nz 1],15);
    f(:,:,:,16)=f(:,[ny 1:ny-1],[nz 1:nz-1],16);
    f(:,:,:,17)=f(:,[ny 1:ny-1],[2:nz 1],17);
    f(:,:,:,18)=f(:,[2:ny 1],[nz 1:nz-1],18);
    f(:,:,:,19)=f(:,[2:ny 1],[2:nz 1],19);
    bouncedback=f(TO_REFLECT); % Densities bouncing back at next timestep
    % Relax; calculate equilibrium state (FEQ) with equivalent speed and density to F
    density = sum(f,4);
    UX=(sum(f(:,:,:,[6 8 9 12 13]),4)-sum(f(:,:,:,[7 10 11 14 15]),4))./density;
    UY=(sum(f(:,:,:,[4 8 10 16 17]),4)-sum(f(:,:,:,[5 9 11 18 19]),4))./density;
    UZ=(sum(f(:,:,:,[2 12 14 16 18]),4)-sum(f(:,:,:,[3 13 15 17 19]),4))./density;
    UX(1,:,:)=UX(1,:,:)+deltaU; %Increase inlet pressure
    UX(on)=0; UY(on)=0; UZ(on)=0; density(on)=0; U_SQU=UX.^2+UY.^2+UZ.^2;
    U8=UX+UY;U9=UX-UY;U10=-UX+UY;U11=-U8;U12=UX+UZ;U13=UX-UZ;
    U14=-U13;U15=-U12;U16=UY+UZ;U17=UY-UZ;U18=-U17;U19=-U16;
    % Calculate equilibrium distribution: stationary
    feq(:,:,:,1) = t1 * density.*(1-3*U_SQU/2);
    % nearest-neighbours
    feq(:,:,:,2) = t2 * density.*(1 + 3*UZ + 9/2*UZ.^2 - 3/2*U_SQU);
    feq(:,:,:,3) = t2 * density.*(1 - 3*UZ + 9/2*UZ.^2 - 3/2*U_SQU);
    feq(:,:,:,4) = t2 * density.*(1 + 3*UY + 9/2*UY.^2 - 3/2*U_SQU);
    feq(:,:,:,5) = t2 * density.*(1 - 3*UY + 9/2*UY.^2 - 3/2*U_SQU);
    feq(:,:,:,6) = t2 * density.*(1 + 3*UX + 9/2*UX.^2 - 3/2*U_SQU);
    feq(:,:,:,7) = t2 * density.*(1 - 3*UX + 9/2*UX.^2 - 3/2*U_SQU);
    % next-nearest neighbours
    feq(:,:,:,8) = t3 * density.*(1 + 3*U8  + 9/2*(U8).^2  - 3*U_SQU/2);
    feq(:,:,:,9) = t3 * density.*(1 + 3*U9  + 9/2*(U9).^2  - 3*U_SQU/2);
    feq(:,:,:,10) = t3 * density.*(1 + 3*U10 + 9/2*(U10).^2 - 3*U_SQU/2);
    feq(:,:,:,11) = t3 * density.*(1 + 3*U11 + 9/2*(U11).^2 - 3*U_SQU/2);
    feq(:,:,:,12) = t3 * density.*(1 + 3*U12 + 9/2*(U12).^2 - 3*U_SQU/2);
    feq(:,:,:,13) = t3 * density.*(1 + 3*U13 + 9/2*(U13).^2 - 3*U_SQU/2);
    feq(:,:,:,14)= t3 * density.*(1 + 3*U14 + 9/2*(U14).^2 - 3*U_SQU/2);
    feq(:,:,:,15)= t3 * density.*(1 + 3*U15 + 9/2*(U15).^2 - 3*U_SQU/2);
    feq(:,:,:,16)= t3 * density.*(1 + 3*U16 + 9/2*(U16).^2 - 3*U_SQU/2);
    feq(:,:,:,17)= t3 * density.*(1 + 3*U17 + 9/2*(U17).^2 - 3*U_SQU/2);
    feq(:,:,:,18)= t3 * density.*(1 + 3*U18 + 9/2*(U18).^2 - 3*U_SQU/2);
    feq(:,:,:,19)= t3 * density.*(1 + 3*U19 + 9/2*(U19).^2 - 3*U_SQU/2);
    f = omega * feq + ( 1 - omega ) * f;
    f(reflected) = bouncedback;
    prevavu=avu;avu=sum(sum(sum(UX)))/numactivenodes; ts=ts+1;
end

%------------plot-------------
figure;
zcut=5;
colormap(hsv(2));
image(2-boundary(:,:,5));
hold on;
quiver(UY(:,:,zcut),UX(:,:,zcut));
xlabel('y');
ylabel('x');
title(['Flow Field at z=',num2str(zcut),', After ',num2str(ts),'\deltat']);

figure;
ycut=5;
colormap(autumn(2));
image(2-squeeze(boundary(:,ycut,:)));
hold on;
quiver(squeeze(UZ(:,ycut,:)),squeeze(UX(:,ycut,:)));
xlabel('z');
ylabel('x');
title(['Flow Field at y=',num2str(ycut),', After ',num2str(ts),'\deltat']);