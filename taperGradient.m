function G = taperGradient(G,lp,l1)
    %% Taper gradient near sources to avoid interference
    [~,n] = size(G);
    taper = sin(linspace(0,pi/2,l1));
    taper = repmat(taper,n,1)';   
    G(lp+1:lp+l1,:)=taper.*G(lp+1:lp+l1,:);
    %% Normalise gradient
    scale = 1.0/max(abs(G(:)));
    G = scale*G; 
end