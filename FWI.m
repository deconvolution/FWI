%% input parameters
clear all;
close all;
load('true_rec.mat');
dims.dh=10;
lpu=100;
lpd=80;
lpl=100;
lpr=100;
v0=ones([51,101])*2000;
v0(15:end,:)=2300;

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
%% Model dimensions
dims.nz=size(v,1);
dims.nx=size(v,2);
% mz=101:151
dims.mz=lpu+1:dims.nz0+lpu;
% mx=101:201
dims.mx=lpl+1:dims.nx0+lpl;
figure('name','initial velocity model');
imagesc(v(dims.mz,dims.mx));
colorbar;
hold on;
plot([lpl+1,dims.nx0+lpl],[lpu+1,lpu+1],'color','red');
hold on;
plot([lpl+1,dims.nx0+lpl],[lpu+dims.nz0,lpu+dims.nz0],'color','red');
hold on;
plot([lpl+1,lpl+1],[lpu+1,lpu+dims.nz0],'color','red');
hold on;
plot([dims.nx0+lpl,dims.nx0+lpl],[lpu+1,lpu+dims.nz0],'color','red');
axis on;
xlabel({['x*' num2str(dims.dh) '[m]']});
ylabel({['z*' num2str(dims.dh) '[m]']});
title('starting velocity model [m/s]');
shg;
%% Receiver locations
dims.rx=min(dims.mx):max(dims.mx);
dims.rz=min(dims.mz)*ones(1,length(dims.rx));
%% Source and source signals
% 21 shot gathers, shootted for every 50m
dims.sz=ones([1,101])*101;
dims.sx=101:201;
% velocity gradient for each source
vg=zeros([size(v),size(dims.sx,2)]);
% times of iteration
n=60;
% source interval during scale searcing
source_interval=1;
n_source_interval=length(1:source_interval:length(dims.sx));

v2=zeros(size(vg,1),size(vg,2),n);
v2(:,:,1)=v;
alp2=zeros([size(vg,3),1]);
of4=zeros([size(alp2,1),1]);
%% top mute
%{
rec_vconst=zeros(dims.nt,length(dims.mx),length(dims.sz));
for i=1:size(dims.sx,2)
    %% forward propagation
    sn=length(dims.sx);
    s0=rickerWave(10,dims);
    s=s0(:,ones(1,length(dims.sz)));
    % some parameters
    nz=dims.nz;
    nx=dims.nx;
    sz=dims.sz(i);
    sx=dims.sx(i);
    dh=dims.dh;
    dt=dims.dt;
    nt=dims.nt;
    rz=dims.rz;
    rx=dims.rx;
    vconst=ones(size(v2(:,:,1)))*v2(1,1,1);
    
    [~,rec_vconst(:,:,i)]=solver2(nz,nx,nt,dh,dt,sz,sx,s(:,i),rz,rx,vconst);
end
%}
%%
tic;
alp=64;
for l=1:n
    fprintf('\n interation=%d/%d',l,n);
    of=zeros([size(dims.sx,2),1]);
    for i=1:size(dims.sx,2)
        %% forward field
        s0=rickerWave(10,dims);
        s=s0(:,ones(1,length(dims.sz)));
        % some parameters
        nz=dims.nz;
        nx=dims.nx;
        sz=dims.sz(i);
        sx=dims.sx(i);
        dh=dims.dh;
        dt=dims.dt;
        nt=dims.nt;
        rz=dims.rz;
        rx=dims.rx;
        
        [p,rec]=solver(nz,nx,nt,dh,dt,sz,sx,s(:,i),rz,rx,v2(:,:,l));
        %fprintf('\nFWI source=%f/%f \t t=%f s',i,size(dims.sz,2),toc);
        %% adjoint field
        % error kernel
        ek=rec-true_rec(:,:,i); %
        % objective function
        of(i)=norm(ek);
        % time reverse error kernel
        ekr=flip(ek,1);
        % inject ekr as the boundary condition back to wave equation
        sz2=dims.rz;
        sx2=dims.rx;
        s2=ekr;
        % p2: adjoint field
        [p2,rec2]=solver(nz,nx,nt,dh,dt,sz2,sx2,s2,rz,rx,v2(:,:,l));
        p2r=flip(p2,3);
        %% calculate gradient
        diff_p2=diff(p2r,1,3)/dt;
        diff_p=diff(p,1,3)/dt;
        for k=1:nt-1
            vg(min(dims.mz):end,:,i)=vg(min(dims.mz):end,:,i)+diff_p2(min(dims.mz):end,:,k).*diff_p(min(dims.mz):end,:,k);
        end
        %fprintf('\nFWI source=%f/%f \t t=%f s',i,size(dims.sz,2),toc);
    end
    of3=mean(of);
    fprintf('\n\t previous OF=%.10f',of3);
    % supress imprint of source
    vg2=taperGradient(sum(vg,3));
    % vg2=sum(vg,3)/(max(max(abs(sum(vg,3)))));
    
    figure(11)
    set(gcf,'Position',[100,100,800,600]);
    subplot(2,2,1)
    imagesc(vg2(dims.mz,dims.mx));
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title('velocity gradient [m/s]');
    axis on;
    shg;
    %% search for scale
    of5=zeros([size(of,1),1]);
    while 1
        fprintf('\n\t alp=%f',alp);
        v2(:,:,l+1)=v2(:,:,l)+alp*vg2;
        for i=1:source_interval:size(dims.sx,2)
            sz=dims.sz(i);
            sx=dims.sx(i);
            [p3,rec3]=solver(nz,nx,nt,dh,dt,sz,sx,s(:,i),rz,rx,v2(:,:,l+1));
            ek3=rec3-true_rec(:,:,i);
            of5(i)=norm(ek3);
            fprintf('\t current OF(%d/%d)=%.10f',i,length(dims.sx),of5(i));
        end
        of6=mean(of5);
        if of6<of3
            alp2(l)=alp;
            break
        else
            alp=alp/2;
            if alp<.001
                break
            end
        end
    end
    if alp<.001
        break
    end
    of4(l)=of6;
    fprintf('\n \t current OF=%.10f',of4(l));
    figure(11)
    subplot(2,2,2)
    imagesc(v2(dims.mz,dims.mx,l+1));
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title('velocity [m/s]');
    axis on;
    subplot(2,2,3)
    scatter(1:l,of4(1:l),'x');
    %xlim([1,l+1]);
    %ylim([min(of4(1:l)),max(of4(1:l))]);
    xlabel('iteration times');
    ylabel('objective function');
    subplot(2,2,4)
    scatter(1:l,alp2(1:l),'x');
    xlabel('iteration times');
    ylabel('alp');
    %xlim([1,l+1]);
    %ylim([alp2(a(1:l)),max(alp2(1:l))]);
    shg;
    %%
    fprintf('\n \t t=%f s',toc);
end
%% plot p
p11=p;
figure(2)
for tt=1:fix(size(p11,3)/20):size(p11,3)
    %C=imfuse(real(p(:,:,i)),v,'falsecolor','Scaling','independent','ColorChannels','red-cyan');
    subplot(2,1,1)
    imshow(p11(:,:,tt),.1*[min(min(p11(:))),max(max(p11(:)))]);
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(tt*dims.dt) 's'],['p']});
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
%{
load('trueModel');
of11=zeros([l,1]);
for i=1:l
    aa=v2(dims.mz,dims.mx,i);
    of11(i)=norm(aa-trueModel);
end
figure('name','objective function');
plot(of11)
%}