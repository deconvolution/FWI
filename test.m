%% input parameters
clear all;
close all;
dims.dh=10;
lpu=100;
lpd=80;
lpl=100;
lpr=100;
load('trueModel.mat');
v0=trueModel;
dims.dt=2*10^-3; % [s]
dims.nz0=size(v0,1); % Cells in z-direction
dims.nx0=size(v0,2); % Cells in x-direction
dims.nt=400; % Amount of time steps
dims.ns=dims.nt;
%% extend the velocity model region
v=v0;
temp=v0(1,:);
temp2=v0(end,:);
v=[temp(ones(lpu,1),:);
    v;
    temp2(ones(lpd,1),:)];
temp3=v(:,1);
temp4=v(:,end);
v=[temp3(:,ones(1,lpl)),v,temp4(:,ones(1,lpr))];
figure;
imagesc(v);
colorbar;
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['z*' num2str(dims.dh) '[m]']});
hold on;
plot([lpl+1,dims.nx0+lpl],[lpu+1,lpu+1],'color','red');
hold on;
plot([lpl+1,dims.nx0+lpl],[lpu+dims.nz0,lpu+dims.nz0],'color','red');
hold on;
plot([lpl+1,lpl+1],[lpu+1,lpu+dims.nz0],'color','red');
hold on;
plot([dims.nx0+lpl,dims.nx0+lpl],[lpu+1,lpu+dims.nz0],'color','red');
axis on;
%% Model dimensions
dims.nz=size(v,1);
dims.nx=size(v,2);
% mz=31:81
dims.mz=lpu+1:dims.nz0+lpu;
% mx=31:131
dims.mx=lpl+1:dims.nx0+lpl;
%% Source and source signals
dims.sz=[100];
dims.sx=[101];
sn=length(dims.sx);
s0=rickerWave(10,dims);
s=s0(:,ones(1,length(dims.sz)));
%% Receiver locations
dims.rx=min(dims.mx):max(dims.mx);
dims.rz=min(dims.mz)*ones(1,length(dims.rx));
%% some parameters
nz=dims.nz;
nx=dims.nx;
sz=dims.sz;
sx=dims.sx;
dh=dims.dh;
dt=dims.dt;
nt=dims.nt;
rz=dims.rz;
rx=dims.rx;

% initial condition
pini=zeros(nz,nx,nt);
%%
[p,rec]=solver(nz,nx,nt,dh,dt,sz,sx,s,rz,rx,v,pini);
%% plot p
figure(2)
for i=1:fix(size(p,3)/20):size(p,3)
    %C=imfuse(real(p(:,:,i)),v,'falsecolor','Scaling','independent','ColorChannels','red-cyan');
    subplot(2,1,1)
    imshow(p(:,:,i),.1*[min(min(p(:))),max(max(p(:)))]);
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(i*dims.dt) 's'],['p']});
    hold on;
    plot([lpl+1,dims.nx0+lpl],[lpu+1,lpu+1],'color','red');
    hold on;
    plot([lpl+1,dims.nx0+lpl],[lpu+dims.nz0,lpu+dims.nz0],'color','red');
    hold on;
    plot([lpl+1,lpl+1],[lpu+1,lpu+dims.nz0],'color','red');
    hold on;
    plot([dims.nx0+lpl,dims.nx0+lpl],[lpu+1,lpu+dims.nz0],'color','red');
    axis on;
    
    subplot(2,1,2)
    imshow(v,[min(v(:)),max(v(:))]);
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    %title({['t=' num2str(i*dims.dt) 's'],['p']});
    hold on;
    plot([lpl+1,dims.nx0+lpl],[lpu+1,lpu+1],'color','red');
    hold on;
    plot([lpl+1,dims.nx0+lpl],[lpu+dims.nz0,lpu+dims.nz0],'color','red');
    hold on;
    plot([lpl+1,lpl+1],[lpu+1,lpu+dims.nz0],'color','red');
    hold on;
    plot([dims.nx0+lpl,dims.nx0+lpl],[lpu+1,lpu+dims.nz0],'color','red');
    axis on;
    shg;
end
%%
figure(3)
imagesc(real(rec),.1*[min(real(rec(:))),max(real(rec(:)))]);
grid on;
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['t*' num2str(dims.dt) '[m]']});
hold off;
%%
t1=real(p(100,:,:));
t2=reshape(t1,[size(t1,2),size(t1,3)]);
for i=1:size(t2,2)
    figure(4)
    subplot(2,1,1)
    plot(t2(:,i));
    ylim([min(t2(:)),max(t2(:))]);
    title(num2str(i));
    shg;
    pause(.01);
end