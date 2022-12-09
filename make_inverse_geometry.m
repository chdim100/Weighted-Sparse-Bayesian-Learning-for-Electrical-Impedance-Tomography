function [imagemodel_pixelized,xc,yc,non_uniform_elems,elementsorted1D,elementsorted2D]...
    =make_inverse_geometry(N,skipcurr,skipvolt,initial_inv_model,I)


[stim, els] = mk_stim_patterns(N,1,[0 skipcurr+1],[0 skipvolt+1],{},I);
initial_inv_model.stimulation=stim;
imagemodel_pixelized=mk_image(initial_inv_model,1);
electrodenodes=[];
for electrode=1:N
    electrodenodes=[electrodenodes initial_inv_model.electrode(electrode).nodes];
end
elements2nodemap=initial_inv_model.elems;
non_uniform_elems=[];
for electode=1:N
    [io,jo]=find(elements2nodemap==electrodenodes(electode));
    non_uniform_elems=[non_uniform_elems; io];
end
non_uniform_elems=unique(non_uniform_elems);
elements2nodemapcleared=elements2nodemap;
elements2nodemapcleared(non_uniform_elems,:)=[];


%%%%%sort triangles to form pixels
total_remained_elems=size(elements2nodemapcleared,1);
remainindex=1:size(elements2nodemapcleared,1);
newindexsorted=[];
nodelist=initial_inv_model.nodes;
ii=1;
%%%%find hypotenuse nodes for each element
hypotenusenodes=[];
for element=1:size(elements2nodemapcleared,1)
    [hypo,hyponodes]=max([norm([nodelist(elements2nodemapcleared(element,1),:)- nodelist(elements2nodemapcleared(element,2),:)])...
        ,norm([nodelist(elements2nodemapcleared(element,2),:)- nodelist(elements2nodemapcleared(element,3),:)])...
        ,norm([nodelist(elements2nodemapcleared(element,1),:)- nodelist(elements2nodemapcleared(element,3),:)])]);
    if hyponodes==1
        hypotenusenodes=[hypotenusenodes; elements2nodemapcleared(element,1) elements2nodemapcleared(element,2)];
    elseif hyponodes==2
        hypotenusenodes=[hypotenusenodes; elements2nodemapcleared(element,2) elements2nodemapcleared(element,3)];
    else
        hypotenusenodes=[hypotenusenodes; elements2nodemapcleared(element,1) elements2nodemapcleared(element,3)];
    end
end

elementsorted1D=[];
elementsorted2D=[];
remainelements=0;
%%%%%find triangle pairs that form each pixel
for element=1:size(elements2nodemapcleared,1)
    if ~ismember(element,elementsorted1D)
        [a, b]=find(hypotenusenodes==hypotenusenodes(element,1));
        [c, d]=find(hypotenusenodes==hypotenusenodes(element,2));
        ind=intersect(a,c);
        ind(ind==element)=[];
        if ~isempty(ind)
            elementsorted2D=[elementsorted2D; element ind];
            elementsorted1D=[elementsorted1D; element; ind];
        end
    end
end

%%%get each pixel's 4 nodes
pixel2nodesmap=zeros(size(elementsorted2D,1),4);
for elem=1:size(elementsorted2D,1)
    pixel2nodesmap(elem,:)=unique([elements2nodemapcleared(elementsorted2D(elem,1),:);...
        elements2nodemapcleared(elementsorted2D(elem,2),:)]);
end

%%%%get pixels' central points
for pixel=1:size(elementsorted2D,1)
    xc(pixel)=mean(nodelist(pixel2nodesmap(pixel,:),1));
    yc(pixel)=mean(nodelist(pixel2nodesmap(pixel,:),2));
end



end