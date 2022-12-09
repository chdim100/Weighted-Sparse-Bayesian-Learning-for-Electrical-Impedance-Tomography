function [dsigma]=iterative_GTR(imagemodel_pixelized,xc,yc,Jpixel,dv,lambda,prior,non_uniform_elems,uniform_elems,elementsorted1D,elementsorted2D,elim,sref,sim)
fprintf('%s:\n', prior)
e=1;
L=size(Jpixel,2);
it=1;   sigmaTi=sref*ones(L,1);
imagemodel_pixelizedk=imagemodel_pixelized;
while e>elim
    if it==1
        Jk=Jpixel;
    else
        [Jfull,Jk]=calculate_pure_Jacobian(imagemodel_pixelizedk,...
            non_uniform_elems,elementsorted1D,1);
    end
    dU=Jk*(sigmaTi-sref);
    %
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
    %
    sigmaTikhnew=(Jk'*Jk+lambda.^2.*Q)\(Jk'*(dv-dU)-lambda.^2.*Q*(sigmaTi-sref));
    imagemodel_pixelized2k=imagemodel_pixelizedk;
    imagemodel_pixelized2k.elem_data(non_uniform_elems)=sref;
    Elementvals=imagemodel_pixelized2k.elem_data;
    Elementvals(non_uniform_elems)=[];
    a=logspace(-4,0,50);
    clear ee
    ee(1)=100;
    for ii=1:50
        Skp1=sigmaTi+a(ii)*sigmaTikhnew;
        %
        if all(Skp1>0)
            ee(ii)=(norm(Jk*(Skp1-sref)-dv))^2+lambda^2*(Skp1-sref)'*Q*(Skp1-sref);
        else
            break;
            %ee(ii)=100;
        end
    end
    if min(ee)>=e
        break;
    else
        e=min(ee);
    end
    fprintf('Iteration %1.0f: cost function value: %1.13f\n',it, min(ee))
    imin=find(ee==e);
    if length(imin)>1
        imin=imin(end);
    end
    sigmaTi=sigmaTi+a(imin)*sigmaTikhnew;
    Elementvals(elementsorted2D(:,1))=sigmaTi;
    Elementvals(elementsorted2D(:,2))=sigmaTi;
    imagemodel_pixelized2k.elem_data(uniform_elems)=Elementvals;
    imagemodel_pixelizedk=imagemodel_pixelized2k;
    it=it+1;
end
dsigma=sigmaTi-sref;
%
% figure
% plot(dv,'b','LineWidth',2)
% hold on
% plot(Jk*dsigma,'r','LineWidth',2)
% legend({'\delta v','J\cdot\delta\sigma'})
% title('\delta v vs J\cdot\delta\sigma')
%

if sim==1
    down=-0.4; up=0.4; step=0.2;
else
    down=-2e-03; up=2e-03; step=1e-03;
end
pp=get_recplot(xc,yc,dsigma,down,up,step);
%
title(['Iterative ' prior])
%
end