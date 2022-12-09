function L=edge_matrix_with_diags(x,y,z)
Lx=length(x);

step=max(diff(x));
e=step/10;
edgesi=1;
if nargin<3
    for pixel=1:Lx
        %x=x(pixel)+step
       right_neighbors{pixel}=find((y==y(pixel)).*(x-x(pixel)<=(step+e)).*(x-x(pixel)>=(step-e)));
       if ~isempty(right_neighbors{pixel})
           edges{edgesi}=[pixel right_neighbors{pixel}];
           edgesi=edgesi+1;
       end
        %y=y(pixel)+step
       upper_neighbors{pixel}=find((x==x(pixel)).*(y-y(pixel)<=(step+e)).*(y-y(pixel)>=(step-e)));
       if ~isempty(upper_neighbors{pixel})
           edges{edgesi}=[pixel upper_neighbors{pixel}];
           edgesi=edgesi+1;
       end
    end
       %upper right diagonal neighbors
       upper_right_neighbors{pixel}=find((y-y(pixel)<=(step+e)).*(y-y(pixel)>=(step-e)).*(x-x(pixel)<=(step+e)).*(x-x(pixel)>=(step-e)));
       if ~isempty(upper_right_neighbors{pixel})
           edges{edgesi}=[pixel upper_right_neighbors{pixel}];
           edgesi=edgesi+1;
       end
       %down right diagonal neighbors
       down_right_neighbors{pixel}=find((y(pixel)-y<=(step+e)).*(y(pixel)-y>=(step-e)).*(x-x(pixel)<=(step+e)).*(x-x(pixel)>=(step-e)));
       if ~isempty(down_right_neighbors{pixel})
           edges{edgesi}=[pixel down_right_neighbors{pixel}];
           edgesi=edgesi+1;
       end
    Nofedges=length(edges);
    L=zeros(Nofedges,Lx);
    for edgesi=1:Nofedges
        nonzeros=edges{edgesi};
        L(edgesi,nonzeros(1))=1;
        L(edgesi,nonzeros(2))=-1;
        if length(nonzeros)>2
        L(edgesi,nonzeros(3))=1;
        L(edgesi,nonzeros(4))=-1;
        end
    end    
else
    error('3D TV is non-available yet')
end
