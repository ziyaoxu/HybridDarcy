function JumpPenaltyEdges = ...
    JumpPenaltyEdges(EtoEmap,BarriersPath,BarriersStartRefCoord,BarriersEndRefCoord)
% time : 2020.6.8 - 2020.6.17
% author : xuziyao
% find out all cell interfaces where we must add penalty on the flux jump (passed by a barrier)
% or we must add penalty on pressure jump (passed by a fracture)
% in other words, 
% find the cell interfaces where we can not add penalty on the pressure jump (passed by a barrier)
% or we can not add penalty on flux jump (passed by a fracture)
GeoTol = 1e-1; % the geometric tolerance that we regard an edge to be aligned with a barrier
JumpPenaltyEdges = zeros(size(EtoEmap));
for m = 1 : length(BarriersPath)
    jj = BarriersPath(m);
    if norm(BarriersStartRefCoord(m,:)-BarriersEndRefCoord(m,:),2) > GeoTol
    if (BarriersStartRefCoord(m,2)+BarriersEndRefCoord(m,2)<GeoTol)&&(EtoEmap(1,jj)~=-1) % aligned with edge 1
        E = EtoEmap(1,jj);
        pp =  EtoEmap(:,E) == jj ;
        JumpPenaltyEdges(1,jj)=1;
        JumpPenaltyEdges(pp,E)=1;
    elseif (2-BarriersStartRefCoord(m,1)-BarriersEndRefCoord(m,1)<GeoTol)&&(EtoEmap(2,jj)~=-1)% aligned with edge 2
        E = EtoEmap(2,jj);
        pp =  EtoEmap(:,E) == jj ;
        JumpPenaltyEdges(2,jj)=1;
        JumpPenaltyEdges(pp,E)=1;
    elseif (2-BarriersStartRefCoord(m,2)-BarriersEndRefCoord(m,2)<GeoTol)&&(EtoEmap(3,jj)~=-1)% aligned with edge 3
        E = EtoEmap(3,jj);
        pp =  EtoEmap(:,E) == jj ;
        JumpPenaltyEdges(3,jj)=1;
        JumpPenaltyEdges(pp,E)=1;
    elseif (BarriersStartRefCoord(m,1)+BarriersEndRefCoord(m,1)<GeoTol)&&(EtoEmap(4,jj)~=-1)% aligned with edge 4
        E = EtoEmap(4,jj);
        pp =  EtoEmap(:,E) == jj ;
        JumpPenaltyEdges(4,jj)=1;
        JumpPenaltyEdges(pp,E)=1;
    else
        for ii = 1 : 4 % the ii-th edge of jj-th element
        E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj
        if E ~= -1 
            pp =  EtoEmap(:,E) == jj ;             
            JumpPenaltyEdges(ii,jj)=1;
            JumpPenaltyEdges(pp,E)=1;
        end
        end
    end
    end
end

end