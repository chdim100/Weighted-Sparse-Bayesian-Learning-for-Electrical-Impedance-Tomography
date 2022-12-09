function w=estimate_weights(xc,yc,g,h,ds,std)
dsr=repmat(abs(ds),[1 h]);
wbar=zeros(g,1);
for i=1:g
    xi=xc(i:i+h-1); yi=yc(i:i+h-1);
    di=get_eucledian_distance2_matrix(xc,yc,xi,yi);
    Gi=get_gaussian(di,std,dsr);
    wbar(i)=sum(sum(Gi));
end
wmax=max(wbar); w=wbar/wmax;
end

function G=get_gaussian(di,std,dsr)
G=dsr.*exp(-di/(2*std^2));
end

function d=get_eucledian_distance2_matrix(xc,yc,xi,yi)
d=(xi-xc').^2+(yi-yc').^2;
end