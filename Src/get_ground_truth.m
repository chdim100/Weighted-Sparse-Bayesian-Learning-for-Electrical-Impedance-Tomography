function reference=get_ground_truth(sim_model,model_code,xc,yc)
%

LC=length(xc);
dist=max(diff(xc));
reference=zeros(LC,1);
LF=length(sim_model.elem_data);
existed=exist('Models\j2d3c\j2d3c.mat');

if strcmp(model_code,'j2d3c')==1&&(length(xc)==1225)&&(existed~=0)
    load('Quick_Data\Circular2Circular\polyelement_inv.mat')
    load('Quick_Data\Circular2Circular\polyelement_sim.mat')
    load('Quick_Data\Circular2Circular\weights.mat')
    for elementinv=1:LC
        for elementsim=1:LF
            reference(elementinv)=reference(elementinv)+w(elementinv,elementsim)*sim_model.elem_data(elementsim);
        end
    end
else
    
    polyelement_inv=zeros(4,2,LC);
    for element=1:LC
        %%%%upper left
        polyelement_inv(1,1,element)=xc(element)-dist/2;
        polyelement_inv(1,2,element)=yc(element)+dist/2;
        %%%%upper right
        polyelement_inv(2,1,element)=xc(element)+dist/2;
        polyelement_inv(2,2,element)=yc(element)+dist/2;
        %%%%down left
        polyelement_inv(4,1,element)=xc(element)-dist/2;
        polyelement_inv(4,2,element)=yc(element)-dist/2;
        %%%%down right
        polyelement_inv(3,1,element)=xc(element)+dist/2;
        polyelement_inv(3,2,element)=yc(element)-dist/2;
    end
    %  
    polyelement_sim=zeros(3,2,LF);
    for element=1:LF
        %%%%%%%
        %%%3 angles per triangle (element)
        for angle=1:3
            polyelement_sim(angle,1,element)=sim_model.fwd_model.nodes(sim_model.fwd_model.elems(element,angle),1);
            polyelement_sim(angle,2,element)=sim_model.fwd_model.nodes(sim_model.fwd_model.elems(element,angle),2);
        end
    end
    
    w=zeros(LC,LF);
    for elementinv=1:LC
        pinverse=polyshape([polyelement_inv(:,1,elementinv) polyelement_inv(:,2,elementinv)]);
        reference(elementinv)=0;
        for elementsim=1:LF
            psim=polyshape([polyelement_sim(:,1,elementsim) polyelement_sim(:,2,elementsim) ]);
            polyout = intersect(pinverse, psim);
            w(elementinv,elementsim)=area(polyout)/area(pinverse);
            reference(elementinv)=reference(elementinv)+w(elementinv,elementsim)*sim_model.elem_data(elementsim);
            if sum(w(elementinv,:))>0.98 %%%speed up if found nearly all corresponding elements!
                %o=[];
                break;
            end
        end
        fprintf('Completed %2.0f /%2.0f\n',elementinv,LC)
    end
    
    if strcmp(model_code,'j2d3c')==1&&(length(xc)==1225)&&(existed==0)
        save('Quick_Data\Circular2Circular\polyelement_inv.mat','polyelement_inv')
        save('Quick_Data\Circular2Circular\polyelement_sim.mat','polyelement_sim')
        save('Quick_Data\Circular2Circular\weights.mat','w')
    end
end
end

