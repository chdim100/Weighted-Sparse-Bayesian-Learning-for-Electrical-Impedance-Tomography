function [coarse_fwd_mdl]=make_2d_geom(invmodel,N,coarseres)

for el=1:N
electrode_nodes(el,:)=invmodel.fwd_model.electrode(el).nodes;
end
electrode_coords=invmodel.fwd_model.nodes(electrode_nodes(:,1),:);
for el=1:N
    elec_nodes{el}=[electrode_coords(el,1) electrode_coords(el,2)];
end



boundarynodes=invmodel.fwd_model.boundary;
boundarycoordinates1=invmodel.fwd_model.nodes(boundarynodes(:,1),:);
boundarycoordinates2=invmodel.fwd_model.nodes(boundarynodes(:,2),:);
b1=1;
for b=1:length(boundarycoordinates1)
    boundarycoordinates(b1,:)=boundarycoordinates1(b,:);
    boundarycoordinates(b1+1,:)=boundarycoordinates2(b,:);
    b1=b1+2;
end

%%%%%make coarse mesh
%coarseres=41;
meshcoarse=min(min(electrode_coords(:,1)),min(electrode_coords(:,2))):2/coarseres:max(max(electrode_coords(:,1)),max(electrode_coords(:,2)));
[X1c,Y1c]=meshgrid(meshcoarse);
xlc=X1c(:);
ylc=Y1c(:);
inc = inpolygon(xlc,ylc,boundarycoordinates(:,1),boundarycoordinates(:,2));
xc=xlc(inc); yc=ylc(inc);

figure
plot(xc,yc,'s','MarkerSize',10.3)
hold on
plot(boundarycoordinates(:,1),boundarycoordinates(:,2),'r')
hold on
plot(electrode_coords(:,1),electrode_coords(:,2),'ko','LineWidth',3)

vtxc= [xc(:),yc(:)];

coarse_fwd_mdl=mk_fmdl_from_nodes(vtxc, elec_nodes,1,'coarse1');
axis square
axis off
% figure
% show_fem(coarse_fwd_mdl)

end

