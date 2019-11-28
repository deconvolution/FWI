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
s_lim=find(abs(source_freq2)<.01*max(abs(source_freq2)) & s_diff<0);
s_lim2=s_lim(1);
f_range=1:s_lim2;
f2=f(f_range);
ome=2*pi*f2;

% source term
sf=source_freq2(1:s_lim2);
%% run solver
nz=dims.nz;
nx=dims.nx;

dh=dims.dh;
nt=dims.nt;
% ome
% sf
rz=dims.rz;
rx=dims.rx;

% v
%% source location in FWI
shot_interval=1;
sx0=dims.sx(1):shot_interval:dims.sx(end);
sz0=dims.sz(1)*ones(size(sx0));
%%

alp=64;
n_iteration=30;
vg_progress=zeros([dims.nz,dims.nx,n_iteration+1]);
v_progress=zeros([dims.nz,dims.nx,n_iteration+1]);
v_progress(:,:,1)=v;
alp_progress=zeros([length(n_iteration)+1,1]);
of_progress=zeros([length(n_iteration)+1,1]);
d0=log(1/Rc)/log(10)*3*v/2/lp;
% wavefield in forward modeling
pt=zeros([dims.nz,dims.nx,dims.nt,length(sx0)]);
% recording in forward modeling
rect=zeros([dims.nt,length(dims.rx),length(sx0)]);
tic;
for l=1:n_iteration
    fprintf('\niteration=%d/%d',l,n_iteration);
    oft=zeros([length(sx0),1]);
    vg_t=zeros([dims.nz,dims.nx,length(sx0)]);
    for i=1:length(sx0)
        % forward modeling
        sz=sz0(i);
        sx=sx0(i);
        if l==1
            [~,pt(:,:,:,i),~,rect(:,:,i)]=solver(nz,nx,dh,nt,sz,sx,ome,n,sf,rz,rx,v_progress(:,:,l),lp,d0);
            et=rect(:,:,i)-true_rect(:,:,dims.sx==sx);
        else
            et=rect(:,:,i)-true_rect(:,:,dims.sx==sx);
        end
        oft(i)=norm(et);
        % time reverse residual
        ert=flip(et,1);
        %% adjoint modeling
        % quantity of time reversed error kernel in frequency domain
        n_er=2000;
        fs2=1/dims.dt;
        f2=fs*(0:(n_er/2))/n_er;
        er=fft(ert,n_er,1)/(n_er/2);
        er2=er(1:n_er/2+1,:);
        t=.01*max(abs(er2),[],1);
        t2=zeros(size(er2,2),1);
        for j=1:length(t2)
            t3=find(abs(er2(:,j))>=t(j));
            t(j)=max(t3);
        end
        f_range2=1:max(t);
        
        sz2=dims.rz;
        sx2=dims.rx;
        
        ome2=2*pi*f2;
        sf2=er2;

        [~,pt2,~,~]=solver(nz,nx,dh,nt,sz2,sx2,ome2(f_range2),n_er,sf2(f_range2,:),rz,rx,v_progress(:,:,l),lp,d0);
        r_pt2=flip(pt2,3);
        r_pt2_t=diff(r_pt2,1,3)/dims.dt;
        pt_t=diff(pt(:,:,:,i),1,3)/dims.dt;
        vg_t(:,:,i)=sum(r_pt2_t.*pt_t,3);
    end
    vg_t2=sum(vg_t,3);
    vg_t3=vg_t2;
    vg_t3(1:dims.mz-1,:)=0;
    vg_progress(:,:,l+1)=taperGradient(vg_t3,lp,l1);
    if l==1
        of_progress(l)=mean(oft);
    end
    fprintf('\n\t previous OF=%.10f',of_progress(l));
    
    figure(11)
    set(gcf,'Position',[100,100,800,600]);
    subplot(2,2,1)
    imagesc(vg_progress(:,:,l+1));
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title('velocity gradient [m/s]');
    axis on;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    hold on;
    ax=plot([lp+1,dims.nx0+lp],[lp+1,lp+1],'color','red');
    hold on;
    ax=plot([lp+1,dims.nx0+lp],[lp+dims.nz0,lp+dims.nz0],'color','red');
    hold on;
    ax=plot([lp+1,lp+1],[lp+1,lp+dims.nz0],'color','red');
    hold on;
    ax=plot([dims.nx0+lp,dims.nx0+lp],[lp+1,lp+dims.nz0],'color','red');
    legend(ax,'PML boundary','location',[0.48,0.03,0.05,0.005],'orientation','horizontal');
    shg;
    %% search for scale
    while 1
        fprintf('\n\t alp=%f',alp);
        v_progress(:,:,l+1)=v_progress(:,:,l)+alp*vg_progress(:,:,l+1);
        d0=log(1/Rc)/log(10)*3*v_progress(:,:,l+1)/2/lp;
        oft3=zeros([length(sx0),1]);
        for i=1:length(sx0)
            sz=sz0(i);
            sx=sx0(i);
            [~,pt(:,:,:,i),~,rect(:,:,i)]=solver(nz,nx,dh,nt,sz,sx,ome,n,sf,rz,rx,v_progress(:,:,l+1),lp,d0);
            et3=rect(:,:,i)-true_rect(:,:,dims.sx==sx);
            oft3(i)=norm(et3);
            fprintf('\t current OF(%d/%d)=%.10f',i,length(sx0),oft3(i));
        end
        of_progress(l+1)=mean(oft3);
        if of_progress(l+1)<of_progress(l)
            alp_progress(l+1)=alp;
            break
        else
            alp=alp/2;
            if alp<1
                break
            end
        end
    end
    if alp<1
        break
    end
    fprintf('\n \t current OF=%.10f',of_progress(l+1));
    figure(11)
    subplot(2,2,2)
    imagesc(v_progress(:,:,l+1));
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title('velocity [m/s]');
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    hold on;
    ax=plot([lp+1,dims.nx0+lp],[lp+1,lp+1],'color','red');
    hold on;
    ax=plot([lp+1,dims.nx0+lp],[lp+dims.nz0,lp+dims.nz0],'color','red');
    hold on;
    ax=plot([lp+1,lp+1],[lp+1,lp+dims.nz0],'color','red');
    hold on;
    ax=plot([dims.nx0+lp,dims.nx0+lp],[lp+1,lp+dims.nz0],'color','red');
    axis on;
    subplot(2,2,3)
    scatter((1:l+1)-1,of_progress(1:l+1),'x');
    xlabel('iteration times');
    ylabel('objective function');
    subplot(2,2,4)
    scatter((1:l+1)-1,alp_progress(1:l+1),'x');
    xlabel('iteration times');
    ylabel('alp');
    shg;
    print(gcf,['C:\Users\zhang\OneDrive\project\acoustic_freq_FWI_PML\progress\progress_image' num2str(l) '.png'],'-dpng','-r600');
    %%
    fprintf('\n \t t=%f s',toc);
end
%% save variable
final.n_iteration=n_iteration;
final.vg=vg_progress;
final.v=v_progress;
final.alp=alp_progress;
final.of=of_progress;
save('final.mat','final');
%% visualize wavefield
%{
ptv=pt;
figure(2)
for i=1:fix(size(pt2,3)/30):size(pt2,3)
    %C=imfuse(real(ptv(:,:,i)),v,'falsecolor','Scaling','independent','ColorChannels','red-cyan');
    subplot(2,1,1)
    imshow(real(ptv(:,:,i)),.1*[min(real(ptv(:))),max(real(ptv(:)))]);
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    title({['t=' num2str(i*dims.dt) 's'],['p']});
    hold on;
    plot([lp+1,dims.nx0+lp],[lp+1,lp+1],'color','red');
    hold on;
    plot([lp+1,dims.nx0+lp],[lp+dims.nz0,lp+dims.nz0],'color','red');
    hold on;
    plot([lp+1,lp+1],[lp+1,lp+dims.nz0],'color','red');
    hold on;
    plot([dims.nx0+lp,dims.nx0+lp],[lp+1,lp+dims.nz0],'color','red');
    axis on;
    
    subplot(2,1,2)
    imshow(vg_progress(:,:,l),[min(min(vg_progress(:,:,l))),max(max(vg_progress(:,:,l)))]);
    colorbar;
    xlabel({['x*' num2str(dims.dh) '[m]']});
    ylabel({['z*' num2str(dims.dh) '[m]']});
    %title({['t=' num2str(i*dims.dt) 's'],['p']});
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
end
%}