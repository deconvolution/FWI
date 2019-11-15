function [p,rec]=solver(nz,nx,nt,dh,dt,sz,sx,s,rz,rx,v)
p=zeros(nz,nx,nt);
s_loc=(sx-1)*nz+sz;
%% index processing
[il,jl]=find(p(:,:,1)==p(:,:,1));
il(s_loc)=[];
jl(s_loc)=[];
dl=find(il==1|il==nz);
dl2=find(jl==1|jl==nx);
dl3=unique([dl;
    dl2]);
il(dl3)=[];
jl(dl3)=[];
%% find non zeros in pini
%{
l_nonzeros=[];
for i=1:size(pini,3)
    if find(pini(:,:,i)~=0)
        l_nonzeros=[l_nonzeros,find(pini(:,:,i)~=0)];
    end
end
l_nonzeros2=unique(l_nonzeros);
%}
%%

% without source
k=1;
for l=1:length(jl)
    p(il(l),jl(l),2)=dt^2*v(il(l),jl(l))^2/2/dh^2*p(il(l),jl(l)+1,1)...
        +dt^2*v(il(l),jl(l))^2/2/dh^2*p(il(l),jl(l)-1,1)...
        +dt^2*v(il(l),jl(l))^2/2/dh^2*p(il(l)+1,jl(l),1)...
        +dt^2*v(il(l),jl(l))^2/2/dh^2*p(il(l)-1,jl(l),1)...
        +(1-2*dt^2*v(il(l),jl(l))^2/dh^2)*p(il(l),jl(l),1);
end
% with source
for ns=1:length(sx)
    p(sz(ns),sx(ns),2)=dt^2*v(sz(ns),sx(ns))^2/2/dh^2*p(sz(ns),sx(ns)+1,1)...
        +dt^2*v(sz(ns),sx(ns))^2/2/dh^2*p(sz(ns),sx(ns)-1,1)...
        +dt^2*v(sz(ns),sx(ns))^2/2/dh^2*p(sz(ns)+1,sx(ns),1)...
        +dt^2*v(sz(ns),sx(ns))^2/2/dh^2*p(sz(ns)-1,sx(ns),1)...
        +(1-2*dt^2*v(sz(ns),sx(ns))^2/dh^2)*p(sz(ns),sx(ns),1)...
        +dt^2/2*s(1,ns);
end
%%
for k=2:nt-1
    % without source
    for l=1:length(jl)
        p(il(l),jl(l),k+1)=dt^2*v(il(l),jl(l))^2/dh^2*p(il(l),jl(l)+1,k)...
            +dt^2*v(il(l),jl(l))^2/dh^2*p(il(l),jl(l)-1,k)...
            +dt^2*v(il(l),jl(l))^2/dh^2*p(il(l)+1,jl(l),k)...
            +dt^2*v(il(l),jl(l))^2/dh^2*p(il(l)-1,jl(l),k)...
            -p(il(l),jl(l),k-1)...
            +(2-4*dt^2*v(il(l),jl(l))^2/dh^2)*p(il(l),jl(l),k);
    end
    % with source
    for ns=1:length(sx)
        p(sz(ns),sx(ns),k+1)=dt^2*v(sz(ns),sx(ns))^2/dh^2*p(sz(ns),sx(ns)+1,k)...
            +dt^2*v(sz(ns),sx(ns))^2/dh^2*p(sz(ns),sx(ns)-1,k)...
            +dt^2*v(sz(ns),sx(ns))^2/dh^2*p(sz(ns)+1,sx(ns),k)...
            +dt^2*v(sz(ns),sx(ns))^2/dh^2*p(sz(ns)-1,sx(ns),k)...
            -p(sz(ns),sx(ns),k-1)...
            +(2-4*dt^2*v(sz(ns),sx(ns))^2/dh^2)*p(sz(ns),sx(ns),k)...
            +dt^2*s(k,ns);
    end
end
%%
rec=reshape(p(rz(1),rx,:),[length(rx),nt]).';
% fprintf('\n solver finished');
end