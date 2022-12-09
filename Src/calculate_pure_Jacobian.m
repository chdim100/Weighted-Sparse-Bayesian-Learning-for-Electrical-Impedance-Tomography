function [Jfull,Jpixel]=...
    calculate_pure_Jacobian(thoracic_imagemodel_pixelized,non_uniform_elems,elementsorted1D,sref)
Jfull=calc_jacobian(thoracic_imagemodel_pixelized);
Jcleared=Jfull;
Jcleared(:,non_uniform_elems)=[];
Jsorted=Jcleared(:,elementsorted1D);
Jpixel=zeros(size(Jsorted,1),size(Jsorted,2)/2);
cj=1;
for coll=1:2:size(Jsorted,2)-1
    Jpixel(:,cj)=(Jsorted(:,coll)+Jsorted(:,coll+1))/2;
    cj=cj+1;
end
Jpixel=Jpixel/sref^2;
%Jnormalized=Jpixel/max(max(abs(Jpixel)));
%Jnormalized=Jpixel;
end