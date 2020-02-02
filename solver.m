function [p,pt,rec,rect]=solver(nz,nx,dh,nt,sz,sx,ome,n,sf,rz,rx,v,lp,d0)
% dims frequency sample number
% Source: source signals in frequency domain
% lp: layers of PML
% d0: for PML, could be related to R (reflection coefficient)
%%
AA=sparse(nz*nx,nz*nx);
s=zeros(nz*nx,length(ome));
p=zeros(nz,nx,length(ome));
%%
for i=1:length(sx)
    s(sz(i)+(sx(i)-1)*nz,:)=sf(:,i);
end
%%
% tic;
%%
for l=2:length(ome)
    %% assign d
    dx=zeros(nz,nx);
    dz=zeros(nz,nx);
    % top
    for i=2:lp+1
        for j=lp+2:nx-lp-1
            dx(i,j)=0;
            dz(i,j)=d0(i,j)*((lp+2-i)/lp)^2;
        end
    end
    % bottom
    for i=nz-lp:nz-1
        for j=lp+2:nx-lp-1
            dx(i,j)=0;
            dz(i,j)=d0(i,j)*((i+1-(nz-lp))/lp)^2;
        end
    end
    % left
    for i=lp+2:nz-lp-1
        for j=2:lp+1
            dz(i,j)=0;
            dx(i,j)=d0(i,j)*((lp+2-j)/lp)^2;
        end
    end
    % right
    for i=lp+1:nz-lp-1
        for j=nx-lp:nx-1
            dz(i,j)=0;
            dx(i,j)=d0(i,j)*((j+1-(nx-lp))/lp)^2;
        end
    end
    % upper left
    for i=2:lp+1
        for j=2:lp+1
            dz(i,j)=d0(i,j)*((lp+2-i)/lp)^2;
            dx(i,j)=d0(i,j)*((lp+2-j)/lp)^2;
        end
    end
    % upper right
    for i=2:lp+1
        for j=nx-lp:nx-1
            dz(i,j)=d0(i,j)*((lp+2-i)/lp)^2;
            dx(i,j)=d0(i,j)*((j+1-(nx-lp))/lp)^2;
        end
    end
    % lower left
    for i=nz-lp:nz-1
        for j=2:lp+1
            dz(i,j)=d0(i,j)*((i+1-(nz-lp))/lp)^2;
            dx(i,j)=d0(i,j)*((lp+2-j)/lp)^2;
        end
    end
    % lower right
    for i=nz-lp:nz-1
        for j=nx-lp:nx-1
            dz(i,j)=d0(i,j)*((i+1-(nz-lp))/lp)^2;
            dx(i,j)=d0(i,j)*((j+1-(nx-lp))/lp)^2;
        end
    end
    gx=1+dx/1i./ome(l);
    gz=1+dz/1i./ome(l);
    gx_x=diff(gx,1,2)/dh;
    gz_z=diff(gz,1,1)/dh;
    gx_x=[-gx_x(:,1:end-31),zeros([size(gx_x,1),1]),gx_x(:,end-30:end)];
    gx_x(:,1)=0;
    gx_x(:,end)=0;
    gz_z=[-gz_z(1:end-31,:);zeros([1,size(gz_z,2)]);gz_z(end-30:end,:)];
    gz_z(1,:)=0;
    gz_z(end,:)=0;
    %% domain interior
    for i=lp+2:nz-1
        for j=lp+2:nx-1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,j
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% top local
    for i=2:lp+1
        for j=lp+2:nx-lp-1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,j
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% left local
    for i=lp+2:nz-lp-1
        for j=2:lp+1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,j
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% upper left
    for i=2:lp+1
        for j=2:lp+1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% upper right
    for i=2:lp+1
        for j=nx-lp:nx-1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,j
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% lower left
    for i=nz-lp:nz-1
        for j=2:lp+1
            % i,j
            AA((j-1)*nz+i,(j-1)*nz+i)=(ome(l)^2-2*v(i,j)^2/gx(i,j)^2/dh^2-2*v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            
            % i,j+1
            AA((j-1)*nz+i,(j)*nz+i)=(v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            % i,j-1
            AA((j-1)*nz+i,(j-2)*nz+i)=(-v(i,j)^2/2/gx(i,j)^3/dh*gx_x(i,j)+v(i,j)^2/gx(i,j)^2/dh^2)/v(i,j)^2;
            
            % i+1,j
            AA((j-1)*nz+i,(j-1)*nz+i+1)=(-v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
            % i-1,j
            AA((j-1)*nz+i,(j-1)*nz+i-1)=(v(i,j)^2/2/gz(i,j)^3/dh*gz_z(i,j)+v(i,j)^2/gz(i,j)^2/dh^2)/v(i,j)^2;
        end
    end
    %% boundary line
    %% B2
    i=1;
    for j=2:nx-1
        AA((j-1)*nz+i,(j-1)*nz+i)=-1/dh-1i*ome(l)/(v(i,j));
        %AA((j-1)*nz+i,(j)*nz+i)=1/dh^2;
        %AA((j-1)*nz+i,(j-2)*nz+i)=1/dh^2;
        %AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh^2;
        AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    end
    %% B3
    i=nz;
    for j=2:nx-1
        AA((j-1)*nz+i,(j-1)*nz+i)=-1/dh-1i*ome(l)/(v(i,j));
        %AA((j-1)*nz+i,(j)*nz+i)=1/dh;
        %AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
        AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
        %AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    end
    %% B4
    j=nx;
    for i=2:nz-1
        AA((j-1)*nz+i,(j-1)*nz+i)=-1/dh-1i*ome(l)/(v(i,j));
        %AA((j-1)*nz+i,(j)*nz+i)=1/dh;
        AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
        %AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
        %AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    end
    %% B5
    j=1;
    for i=2:nz-1
        AA((j-1)*nz+i,(j-1)*nz+i)=-1/dh-1i*ome(l)/(v(i,j));
        AA((j-1)*nz+i,(j)*nz+i)=1/dh;
        %AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
        %AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
        %AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    end
    %% B6
    j=nx;
    i=1;
    AA((j-1)*nz+i,(j-1)*nz+i)=-2/dh-2i*ome(l)/(v(i,j));
    %AA((j-1)*nz+i,(j)*nz+i)=1/dh;
    AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
    %AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
    AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    %% B7
    j=1;
    i=1;
    AA((j-1)*nz+i,(j-1)*nz+i)=-2/dh-2i*ome(l)/(v(i,j));
    AA((j-1)*nz+i,(j)*nz+i)=1/dh;
    %AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
    %AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
    AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    %% B8
    i=nz;
    j=nx;
    AA((j-1)*nz+i,(j-1)*nz+i)=-2/dh-2i*ome(l)/(v(i,j));
    %AA((j-1)*nz+i,(j)*nz+i)=1/dh;
    AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
    AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
    %AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    %% B9
    j=1;
    i=nz;
    AA((j-1)*nz+i,(j-1)*nz+i)=-2/dh-2i*ome(l)/(v(i,j));
    AA((j-1)*nz+i,(j)*nz+i)=1/dh;
    %AA((j-1)*nz+i,(j-2)*nz+i)=1/dh;
    AA((j-1)*nz+i,(j-1)*nz+i-1)=1/dh;
    %AA((j-1)*nz+i,(j-1)*nz+i+1)=1/dh;
    %% solve Ax=b
    p0=AA\(-s(:,l));
    p(:,:,l)=reshape(p0,[nz,nx]);
    
    %fprintf('frequency step=%f/%f \t',l,length(ome));
    %fprintf('t=%f s\n',toc);
end
%%
rec=reshape(p(rz(1),rx,:),[length(rx),length(ome)]).';
pt0=ifft(p,n,3)*n;
pt=real(pt0(:,:,1:nt));
rect=reshape(pt(rz(1),rx,:),[length(rx),nt]).';
end