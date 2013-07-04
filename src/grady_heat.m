%% Setup surface, nodes, and radial kernel
% Red blood cell (t2) surface
r0=1;
r1=2;
t2 = @(la,th) [(r1 + r0*cos(th)).*cos(la) (r1+r0*cos(th)).*sin(la) r0*sin(th)];
syms la th; % Compute the surface normal vectors symbolically
nr = inline(cross(diff(t2(la,th),la),diff(t2(la,th),th)));

% Use minimum energy nodes projected to the sphere
x=getMEPoints(20);
x=x(:,1:3);
N=size(x,1);

[la,th]=cart2sph(x(:,1),x(:,2),x(:,3));
x=t2(la,th);
nr=nr(la,th); nr=nr./repmat(sqrt(sum(nr.^2,2)),[1 3]);
% Radial kernel
ep=4;  % Shape parameter
phi=@(r2) 1./sqrt(1+ep^2*r2);         % IMQ
dphi=@(r2) -ep^2./sqrt(1+ep^2*r2).^3; % Derivative of IMQ over r

%%  Compute the Surface Laplacian - not specific to t2 surface
xij=repmat(x(:,1),[1 N]); xij=xij-xij.'; nxi=repmat(nr(:,1),[1 N]);
yij=repmat(x(:,2),[1 N]); yij=yij-yij.'; nyi=repmat(nr(:,2),[1 N]);
zij=repmat(x(:,3),[1 N]); zij=zij-zij.'; nzi=repmat(nr(:,3),[1 N]);
r2=xij.^2 + yij.^2 + zij.^2; A=dphi(r2);
DPx=((1-nxi.^2).*xij - nxi.*nyi.*yij - nxi.*nzi.*zij).*A;
DPy=(-nxi.*nyi.*xij + (1-nyi.^2).*yij - nyi.*nzi.*zij).*A;
DPz=(-nxi.*nzi.*xij - nyi.*nzi.*yij + (1-nzi.^2).*zij).*A;
A=chol(phi(r2)); DPx=(DPx/A)/A.'; DPy=(DPy/A)/A.'; DPz=(DPz/A)/A.';
Lap=DPx*DPx+DPy*DPy+DPz*DPz;  % Surface Laplacian

%% Turing spot pattern example using SBDF2
del=0.0045;
d=0.516;
tau1=0.02;
tau2=0.2;
alp=0.899;
bet=-0.91;
gam=-alp;
ufun=@(u,v) alp*u.*(1-tau1*v.^2)+v.*(1-tau2*u);
vfun=@(u,v) bet*v.*(1+alp*tau1/bet*u.*v)+u.*(gam+tau2*v);

% Initial condition
u=u0;


% Implicit systems for SBDF2
dt=0.05;
tfinal=200;
Du=1.5*eye(N)-del*d*dt*(Lap);
[Lu,Uu,pu]=lu(Du,'vector');
Dv=1.5*eye(N)-del*dt*(Lap);
[Lv,Uv,pv]=lu(Dv,'vector');
optsu.UT=true;
optsl.LT=true;
% One step of backward Euler to bootstrap SBDF2
rhsu=u;
u0=u;
u=(eye(N)-del*dt*Lap)\rhsu;

sz=[101 201];
M = prod(sz);  % surface grid parameters
[ll,tt]=meshgrid(linspace(-pi,pi,sz(2)),linspace(-pi,pi,sz(1)));
xx=t2(ll(:),tt(:));
% Interpolate to the grid
re2=(repmat(xx(:,1),[1 N])-repmat(x(:,1).',[M 1])).^2;
re2=re2+(repmat(xx(:,2),[1 N])-repmat(x(:,2).',[M 1])).^2;
re2=re2+(repmat(xx(:,3),[1 N])-repmat(x(:,3).',[M 1])).^2;

yy=reshape(xx(:,2),sz);
zz=reshape(xx(:,3),sz);
xx=reshape(xx(:,1),sz);
    
for j=1:tfinal/dt
    rhsu=u;
    u0=u;
    u=(eye(N)-del*dt*Lap)\rhsu;
    
    
    uu=reshape(phi(re2)*(A\(A.'\u)),sz); 
    % Plot the results
    h=surf(xx,yy,zz,uu);
    shading interp;
    daspect([1 1 1]);
    axis tight;
    
    %caxis([gmin gmax]);
    set(h, 'edgecolor','none')
    set(gca,'Visible','off')
    set(gcf,'PaperPositionMode','auto')

    name = sprintf('../output4/imageUNSCALED_%03d.png',j-1);
    print('-dpng',name);
    close all force;
    close all hidden;
    
end

