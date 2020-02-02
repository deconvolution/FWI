%% input parameters
clear all;
close all;
load('true_rect.mat');
dims.dh=10;
lp=30;
% initial velocity model
v0=ones([51,101])*2000;
% estimation to first layer thickness
l1=15;
% estimation to first layer velocity
v0(l1:end,:)=2300;
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
figure;
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
%% Model dimensions
dims.nz=size(v,1);
dims.nx=size(v,2);
% mz=31:81
dims.mz=lp+1:dims.nz0+lp;
% mx=31:131
dims.mx=lp+1:dims.nx0+lp;
%% source
singles=rickerWave(10,dims);
dims.sx=min(dims.mx):max(dims.mx);
dims.sz=min(dims.mz)*ones(1,length(dims.mx));
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
title('frequency domain');
shg;
%% find effective frequency of source
s_diff=diff(abs(source_freq2));
s_diff=[s_diff(1);s_diff];
s_lim=find(abs(source_freq2)<.1*max(abs(source_freq2)) & s_diff<0);
s_lim2=s_lim(1);
f_range=1:s_lim2;
f2=f(f_range);
ome=2*pi*f2;
% source term
sf=source_freq2(1:s_lim2);
%% source location in FWI
shot_interval=10;
sx0=dims.sx(1):shot_interval:dims.sx(end);
sz0=dims.sz(1)*ones(size(sx0));
% extend source term
sf2=sf*ones(1,length(sx0));
% sf2=sf(:,ones(3,1));
%% find effective frequency for true_rect
n_true_recf=n;
fs2=fs;
f2=fs2*(0:(n_true_recf/2))/n_true_recf;
true_recf=fft(true_rect,n_true_recf,1)/(n_true_recf/2);
true_recf=true_recf(1:s_lim2,:,:);
%% solver parameters
nz=dims.nz;
nx=dims.nx;

dh=dims.dh;
nt=dims.nt;
% ome
% sf
rz=dims.rz;
rx=dims.rx;

n_iteration=20;
alp_progress=zeros([length(n_iteration)+1,1]);
of_progress=zeros([length(n_iteration)+1,1]);

d0=log(1/Rc)/log(10)*3*v/2/lp;
gmax1=150;
C=zeros(length(n_iteration),1);
FID=fopen('log.txt','w');
tic;
for l=1:n_iteration
    fprintf('\n current iteration=%d/%d \t time=%f s',l,n_iteration,toc);
    fprintf(FID,'\n current iteration=%d/%d \t time=%f s',l,n_iteration,toc);
    %%
    vg=zeros(size(v));
    %% gradient
    for i=1:length(sx0)
        loc=(i-1)*shot_interval+1;
        fprintf('\n \t current shot number=%d/%d',i,length(sx0));
        fprintf(FID,'\n \t current shot number=%d/%d',i,length(sx0));
        true_recfs=true_recf(:,:,loc);
        [pf,recf,Ct,vgt]=FWI_solver(nz,nx,dh,sz0(i),sx0(i),ome,sf2(:,i),rz,rx,v,lp,d0,true_recfs);
        vg=vgt+vg;
    end
    vg=vg/length(sx0);
    vg=taperGradient(vg,lp,l1);
    C(l)=Ct/length(sx0);
    fprintf('\n \t\t\t Cost=%f',C(l));
    fprintf(FID,'\n \t Cost=%f',C(l));
    fprintf(FID,'\n \t max gradient=%f',max(abs(vg(:))));
    %% optimization
    % max gradient for the first iteration
    if l==1
        alp=gmax1/max(abs(vg(:)));
    end
    v=v-alp*vg;
    figure(99)
    subplot(2,1,1)
    imagesc(vg);
    hold on;
    plot([lp+1,nx-lp],[lp+1,lp+1],'color','red');
    hold on;
    plot([lp+1,nx-lp],[nz-lp,nz-lp],'color','red');
    hold on;
    plot([lp+1,lp+1],[lp+1,nz-lp],'color','red');
    hold on;
    plot([nx-lp,nx-lp],[lp+1,nz-lp],'color','red');
    title({['G [m/s]'],[ 'max gradient= ' num2str(max(abs(vg(:))))]});
    colorbar;
    subplot(2,1,2)
    imagesc(v);
    hold on;
    plot([lp+1,nx-lp],[lp+1,lp+1],'color','red');
    hold on;
    plot([lp+1,nx-lp],[nz-lp,nz-lp],'color','red');
    hold on;
    plot([lp+1,lp+1],[lp+1,nz-lp],'color','red');
    hold on;
    plot([nx-lp,nx-lp],[lp+1,nz-lp],'color','red');
    colorbar;
    title('v [m/s]');
    shg;
end
fclose(FID);
%%
figure('name','final result');
subplot(2,2,1)
load('trueModel.mat');
imagesc(trueModel);
colorbar;
title('true model');
subplot(2,2,2)
imagesc(v(dims.mz,dims.mx));
colorbar;
title('inversion result');
subplot(2,2,3)
imagesc(v(dims.mz,dims.mx)-trueModel);
title('error');
%%
