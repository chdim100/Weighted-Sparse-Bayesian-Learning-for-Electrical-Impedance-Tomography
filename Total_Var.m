function [dsigmaTV]=Total_Var(imagemodel_pixelized,xc,yc,Jpixel,dv,lambda,bhta,cofactor,non_uniform_elems,uniform_elems,elementsorted1D,elementsorted2D,elim,maxIT,sref,sim)
fprintf('Total Variation:\n')
e=1; 
L=size(Jpixel,2);
it=1;   sigmaTV=sref*ones(L,1);
imagemodel_pixelizedk=imagemodel_pixelized;
Led=edge_matrix_with_diags(xc,yc);
nofedges=size(Led,1);
xu=zeros(nofedges,1);
while e>elim&&it<=maxIT
    if it==1
        Jk=Jpixel;
    else
        [Jfull,Jk]=calculate_pure_Jacobian(imagemodel_pixelizedk,...
        non_uniform_elems,elementsorted1D,1);
    end
    dsigmaTV=sigmaTV-sref;
    dU=Jk*dsigmaTV;
    
    %%%%%%%
    bhta=bhta/cofactor;
    Ek=zeros(nofedges,nofedges);
    K=zeros(nofedges,nofedges);
    for edgei=1:nofedges
        Ek(edgei,edgei)=sqrt((norm(Led(edgei,:)*dsigmaTV))^2+bhta);
        K(edgei,edgei)=1-xu(edgei)*Led(edgei,:)*dsigmaTV/Ek(edgei,edgei);
    end
    Einv=diag(diag(Ek).^-1);
    %%%% updates of dual variables
    sigmaTVkhnew=left_divide((Jk'*Jk+lambda^2*Led'*Einv*K*Led),...
        (Jk'*(dv-dU)-lambda^2*Led'*Einv*Led*dsigmaTV));
    xu=-xu+Einv*Led*dsigmaTV+Einv*K*Led*sigmaTVkhnew;
    xu=dual_variable_scaling_rule(xu);
    %%%%%%%

    imagemodel_pixelized2k=imagemodel_pixelizedk;
    imagemodel_pixelized2k.elem_data(non_uniform_elems)=sref;
    Elementvals=imagemodel_pixelized2k.elem_data;
    Elementvals(non_uniform_elems)=[];
    a=logspace(-4,0,50);
    for ii=1:50
        Skp1=sigmaTV+a(ii)*sigmaTVkhnew;
        TVterm=0;
        for edgei=1:nofedges
            TVterm=TVterm+sqrt((norm(Led(edgei,:)*(Skp1-sref)))^2+bhta);
        end
        if all(Skp1>0)
           ee(ii)=0.5*(norm(Jk*(Skp1-sref)-dv))^2+0.5*lambda^2*TVterm;
        else
            break;
            %ee(ii)=100;
        end
            end
    fprintf('TV cost function value (Iteration %1.0f): %1.13f\n', it, min(ee))
    %
    if min(ee)>=e
        break;
    else
    e=min(ee);
    end
    imin=find(ee==e);
    if length(imin)>1
        imin=imin(end);
    end
    sigmaTV=sigmaTV+a(imin)*sigmaTVkhnew;
    Elementvals(elementsorted2D(:,1))=sigmaTV;
    Elementvals(elementsorted2D(:,2))=sigmaTV;
    imagemodel_pixelized2k.elem_data(uniform_elems)=Elementvals;
    imagemodel_pixelizedk=imagemodel_pixelized2k;
    it=it+1;
end
%
dsigmaTV=sigmaTV-sref;
% figure
% plot(dv,'b','LineWidth',2)
% hold on
% plot(Jk*dsigmaTV,'r','LineWidth',2)
% legend({'\delta v','J\cdot\delta\sigma'})
% title('\delta v vs J\cdot\delta\sigma')
%
%
if sim==1
    down=-0.4; up=0.4; step=0.2;
else
    down=-2e-03; up=2e-03; step=1e-03;
end
pp=get_recplot(xc,yc,dsigmaTV,down,up,step);
title('TV Iterative')

end