function dsigma=single_step_GTR(Jpixel,xc,yc,dv,prior,lambda,sim)
%
L=size(Jpixel,2);
switch prior
    case {'STD','Standard Tikhonov', 'Tikhonov'}
        Q=eye(L,L);
    case {'NOSER','Noser','noser'}
        Q=diag(diag(Jpixel'*Jpixel));
    case {'Laplace','laplace'}
        Q=calc_pix_Laplacian(xc,yc);
    otherwise
        error('Not Included')
end
dsigma=(Jpixel'*Jpixel+lambda^2*Q)\(Jpixel'*dv);
%
if sim==1
    down=-0.4; up=0.4; step=0.2;
else
    down=-0.003; up=0.003; step=0.001;
end
pp=get_recplot(xc,yc,dsigma,down,up,step);
title(['Single-Step ' prior])
%
end