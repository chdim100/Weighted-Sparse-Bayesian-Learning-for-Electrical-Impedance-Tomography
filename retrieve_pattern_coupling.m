function Correlations=retrieve_pattern_coupling(xc,yc,h)
    g=length(xc)-h+1;
    for i=1:g
        [iprev,inext]=find_neighbour_clusters(i,h,xc,yc);
        Correlations(i,1)=iprev;
        Correlations(i,2)=inext;
    end
end