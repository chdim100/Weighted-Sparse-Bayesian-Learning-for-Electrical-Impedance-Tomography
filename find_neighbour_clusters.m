function [iprev,inext]=find_neighbour_clusters(i,h,xc,yc)

%i=43; h=4;

step=xc(11)-xc(10); ee=step/10; iprev=-1; inext=-1;
%%%%find neigbhour clusters to perform pattern coupling
xi=xc(i:i+h-1); yi=yc(i:i+h-1);

for ipoint=1:h
    indices_of_prev{ipoint}=find(xc>=xi(ipoint)-step-ee&xc<=xi(ipoint)-step+ee&yc==yi(ipoint));
    if isempty(indices_of_prev{ipoint})
        iprev=-2;
        break;
    end
end
%
if iprev~=-2
    pointsprev=[indices_of_prev{1:h}];
    if ismember(pointsprev,indices_of_prev{1}:indices_of_prev{1}+h-1)
        iprev=indices_of_prev{1};
    end
end

for ipoint=1:h
    indices_of_next{ipoint}=find(xc>=xi(ipoint)+step-ee&xc<=xi(ipoint)+step+ee&yc==yi(ipoint));
    if isempty(indices_of_next{ipoint})
        inext=-2;
        break;
    end
end


if inext~=-2
    pointsnext=[indices_of_next{1:h}];
    if ismember(pointsnext,indices_of_next{1}:indices_of_next{1}+h-1)
        inext=indices_of_next{1};
    end
end

% %%%%%
% figure
% plot(xc,yc,'b*')
% hold on
% plot(xi,yi,'ro')
% hold on
% plot(xc(indices_of_next{1}),yc(indices_of_next{1}),'ks')
% hold on
% plot(xc(indices_of_next{2}),yc(indices_of_next{2}),'ks')
% hold on
% plot(xc(indices_of_next{3}),yc(indices_of_next{3}),'ks')
% hold on
% plot(xc(indices_of_next{4}),yc(indices_of_next{4}),'ks')
% % %%%%%

end