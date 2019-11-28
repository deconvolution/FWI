%% input parameters
clear all;
close all;
dims.dh=10;
lp=30;
load('trueModel.mat');
v0=trueModel;
dims.dt=10^-3; % [s]
dims.nz0=size(v0,1); % Cells in z-direction
dims.nx0=size(v0,2); % Cells in x-direction
dims.nt=1000; % Amount of time steps
dims.ns=dims.nt;
%% extend the velocity model region to employ PML
v=v0;
temp=v0(1,:);
temp2=v0(end,:);
v=[temp(ones(lp,1),:);
    v;
    temp2(ones(lp,1),:)];
temp3=v(:,1);
temp4=v(:,end);
v=[temp3(:,ones(1,lp)),v,temp4(:,ones(1,lp))];
figure('name','true velocity');
imagesc(v);
colorbar;
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['z*' num2str(dims.dh) '[m]']});
hold on;
plot([lp+1,dims.nx0+lp],[lp+1,lp+1],'color','red');
hold on;
plot([lp+1,dims.nx0+lp],[lp+dims.nz0,lp+dims.nz0],'color','red');
hold on;
plot([lp+1,lp+1],[lp+1,lp+dims.nz0],'color','red');
hold on;
plot([dims.nx0+lp,dims.nx0+lp],[lp+1,lp+dims.nz0],'color','red');
axis on;
shg;
Rc=.1;
d0=log(1/Rc)/log(10)*3*v/2/lp;
%% Model dimensions
dims.nz=size(v,1);
dims.nx=size(v,2);
% mz=31:81
dims.mz=lp+1:dims.nz0+lp;
% mx=31:131
dims.mx=lp+1:dims.nx0+lp;
%% source signals
singles=rickerWave(10,dims);
%% Receiver locations
dims.rx=min(dims.mx):max(dims.mx);
dims.rz=min(dims.mz)*ones(1,length(dims.rx));
%% source
fs=1/dims.dt;
L=dims.nt;
n=2000;
f=fs*(0:(n/2))/n;
s=zeros([dims.nz*dims.nx,1]);
source_freq=fft(singles,n)/(n/2);
source_freq2=source_freq(1:n/2+1);
% check source
source_time=ifft(source_freq2,n,1)*n;
figure('name','source signals');
subplot(3,1,1)
plot(dims.dt:dims.dt:dims.dt*L,real(singles));
xlabel(['t [s]']);
title('original');
subplot(3,1,2)
plot(dims.dt:dims.dt:dims.dt*L,real(source_time(1:dims.nt)));
xlabel(['t' [s]']);
title('transformed from frequency');
subplot(3,1,3)
plot(abs(source_freq2));
%aa=reshape(aa,[300,1]);
%plot(real(aa))
title('frequencu domain');
shg;
%% find effective frequency of source
s_diff=diff(abs(source_freq2));
s_diff=[s_diff(1);s_diff];
s_lim=find(abs(source_freq2)<.01*max(abs(source_freq2)) & s_diff<0);
s_lim2=s_lim(1);
f_range=1:s_lim2;
f2=f(f_range);
ome=2*pi*f2;

% source term
sf=source_freq2(1:s_lim2);

source_time2=ifft(source_freq2(f_range),n,1)*n;
figure('name','check source')
plot(dims.dt:dims.dt:dims.dt*L,real(source_time2(1:dims.nt)));
%% run solver
nz=dims.nz;
nx=dims.nx;

dh=dims.dh;
nt=dims.nt;
% ome
% sf
rz=dims.rz;
rx=dims.rx;

% boundary condition
pini=zeros(nz,nx,length(ome));

% v
%% source location
sx0=dims.mx(1):1:dims.mx(end);
sz0=dims.mz(1)*ones(size(sx0));

% true_rec: (time_samples,receivers,source), common shot gather
true_rect=zeros(dims.nt,length(rx),length(sx0));
%%
tic;
for i=1:length(sx0)
    sz=sz0(i);
    sx=sx0(i);
    [~,pt,~,rect]=solver(nz,nx,dh,nt,sz,sx,ome,n,sf,rz,rx,v,lp,d0);
    true_rect(:,:,i)=rect;
    fprintf('\ncommon shot gather generator \t source=%f/%f \t',i,length(sx0));
    fprintf('t=%f s', toc);
end
%%
save('true_rect.mat','true_rect');
%%
load('true_rect.mat');
for i=1:size(true_rect,3)
    figure(5)
    imagesc(true_rect(:,:,i),1*[min(true_rect(:)),max(true_rect(:))]);
    shg;
    pause(.5);
end