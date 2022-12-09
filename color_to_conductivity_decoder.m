function [ImCond,xflatten,yflatten]=color_to_conductivity_decoder(perturbations,datapath)
sel1=perturbations(1);
sel2=perturbations(2);
filename=[datapath '\Edited_photos\fantom_',num2str(sel1),'_',num2str(sel2),'_simple.jpg'];
%[Image,~]=imread(filename,'jpg');
[Image,~]=imread(filename);
Mindist=min(size(Image,1),size(Image,2));
Image=Image(1:Mindist,1:Mindist,:);
xim=linspace(-1,1,Mindist);
yim=(linspace(1,-1,Mindist));
threshold=30;
ImCond=ones(size(Image,1),size(Image,2));
for i=1:size(Image,1)
    for j=1:size(Image,2)
        if Image(i,j,1)>=255-threshold&&Image(i,j,2)<threshold&&Image(i,j,3)<threshold
            ImCond(i,j)=10;
        elseif Image(i,j,1)<57+threshold&&Image(i,j,2)>=255-threshold&&Image(i,j,3)<threshold
            ImCond(i,j)=1;
        elseif Image(i,j,1)<threshold&&Image(i,j,2)<threshold&&Image(i,j,3)>=255-threshold
            ImCond(i,j)=0.1;
        elseif xim(i)^2+yim(j)^2>1
            ImCond(i,j)=-1;
        end
    end
end
[xflatten,yflatten]=meshgrid(xim,yim);
xflatten=xflatten(:);
yflatten=yflatten(:);
ImCond=ImCond(:);
indices=find(ImCond~=-1);
xflatten=xflatten(indices);
yflatten=yflatten(indices);
ImCond=ImCond(indices);
end