function [coordinates,elements,EtoEmap] = ElementMesh( xmin,xmax,ymin,ymax,N,M )
%   本函数目前仅适用于矩形区域的均匀等腰直角三角形剖分.
XX=linspace(xmin,xmax,N+1); YY=linspace(ymin,ymax,M+1); 
[X,Y]=meshgrid(XX,YY);      % 纵线横线交错生成方格网节点的x、y坐标，分别存入X、Y矩阵
coordinates=zeros((N+1)*(M+1),2); elements=zeros(3,2*N*M);  
X=X'; Y=Y';    
coordinates(:,1)=X(:);coordinates(:,2)=Y(:); 
for I=1:M            
    for J=1:N
        star=(I-1)*(N+1)+J;    
        num=((I-1)*N+J-1)*2+1; 
        elements(:,num)=[star,star+1,star+N+1]; 
        elements(:,num+1)=[star+1,star+N+2,star+N+1];
    end 
end 
EtoEmap = triangulation_neighbor_triangles ( elements ); 

Geoeps = ( abs(xmax-xmin)/N + abs(ymax-ymin)/M ) * 1e-6 ; % Geometric tolerance
boundary_edge = triangulation_order3_boundary_edge ( elements );
boundary_edge_midp = ( coordinates(boundary_edge(1,:)',:) + ...
    coordinates(boundary_edge(2,:)',:) ) / 2 ;
another_edge = zeros( size(boundary_edge_midp,1),1 );
for i = 1 : size(boundary_edge_midp,1)
    xi = boundary_edge_midp(i,1);
    yi = boundary_edge_midp(i,2);
    for j=1 : size(boundary_edge_midp,1)
    xj = boundary_edge_midp(j,1);
    yj = boundary_edge_midp(j,2);
        if ((abs(abs(xi-xj)-abs(xmax-xmin))<Geoeps)&&(abs(yi-yj)<Geoeps))...
                || ((abs(abs(yi-yj)-abs(ymax-ymin))<Geoeps)&&(abs(xi-xj)<Geoeps))
            another_edge(i)=j;
            continue;
        end 
    end
end
edge_to_element = zeros(size(boundary_edge,2),1);
boundary_element = find(min(EtoEmap)<0)';
for i = 1 : size(edge_to_element)
    for j = 1: size(boundary_element)
        if all(ismember(boundary_edge(:,i),elements(:,boundary_element(j))))
            edge_to_element(i) = boundary_element(j) ; 
            continue ; 
        end
    end
end
for j = 1 : size(boundary_element)
   e = elements(:,boundary_element(j))';
    if EtoEmap(1,boundary_element(j)) == (-1)
        EtoEmap(1,boundary_element(j))=...
            edge_to_element(another_edge(...
            all(ismember(boundary_edge',e([2,3])),2)>0));
    end
    if EtoEmap(2,boundary_element(j)) == (-1)
        EtoEmap(2,boundary_element(j))=...
            edge_to_element(another_edge(...
            all(ismember(boundary_edge',e([3,1])),2)>0));
    end
    if EtoEmap(3,boundary_element(j)) == (-1)
        EtoEmap(3,boundary_element(j))=...
            edge_to_element(another_edge(...
            all(ismember(boundary_edge',e([1,2])),2)>0));
    end
end
end