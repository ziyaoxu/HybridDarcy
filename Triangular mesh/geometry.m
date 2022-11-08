Xmin   = 0   ; % left boundary of the Domain 
Xmax   = 10  ; % right boundary of the Domain 
Ymin   = 0   ; % lower boundary of the Domain 
Ymax   = 10  ; % upper boundary of the Domain 
Geoeps = ( abs(Xmax-Xmin) + abs(Ymax-Ymin) ) * 1e-7 ; % Geometric tolerance
h0=2 ; % Initial edge length 
xvec = linspace(Xmin,Xmax,int32(1+(Xmax-Xmin)/h0))';
yvec = linspace(Ymin,Ymax,int32(1+(Ymax-Ymin)/h0))';
xvec(1)=[];xvec(end)=[];yvec(1)=[];yvec(end)=[];
pfix=[[Xmin,Ymin];... % boundary node is Fixed (NFIXx2)
    [xvec,Ymin*ones(size(xvec))];...
    [Xmax,Ymin];...
    [Xmax*ones(size(yvec)),yvec];...
    [Xmax,Ymax];...
    [xvec(end:-1:1),Ymax*ones(size(xvec))];...
    [Xmin,Ymax];...
    [Xmin*ones(size(yvec)),yvec(end:-1:1)]];
fd=@(p) drectangle(p,Xmin,Xmax,Ymin,Ymax); 
[coordinates,t10086]=distmesh2d(fd,@huniform,h0,[Xmin,Ymin;Xmax,Ymax],pfix);
elements = transpose(t10086) ;
boundary_edge = triangulation_order3_boundary_edge ( elements );
EtoEmap = triangulation_neighbor_triangles ( elements ); 
boundary_edge_midp = ( coordinates(boundary_edge(1,:)',:) + ...
    coordinates(boundary_edge(2,:)',:) ) / 2 ;
another_edge = zeros( size(boundary_edge_midp,1),1 );
for i = 1 : size(boundary_edge_midp,1)
    xi = boundary_edge_midp(i,1);
    yi = boundary_edge_midp(i,2);
    for j=1 : size(boundary_edge_midp,1)
    xj = boundary_edge_midp(j,1);
    yj = boundary_edge_midp(j,2);
        if ((abs(abs(xi-xj)-abs(Xmax-Xmin))<Geoeps)&&(abs(yi-yj)<Geoeps))...
                || ((abs(abs(yi-yj)-abs(Ymax-Ymin))<Geoeps)&&(abs(xi-xj)<Geoeps))
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
%clearvars -EXCEPT  coordinates elements EtoEmap 
%close all
triangle_node = elements' ;
% 写入节点信息：
  output_unit = fopen ( 'MG_nodes.txt', 'wt' );
  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'Error!' );
  end
  for ii = 1 : size(coordinates,1)
    for jj = 1 : size(coordinates,2)
      fprintf ( output_unit, '  %14f', coordinates(ii,jj) );
    end
    fprintf ( output_unit, '\n' );
  end
  fclose ( output_unit );
% 写入三角元信息：
  output_unit = fopen ( 'MG_elements.txt', 'wt' );
  if ( output_unit < 0 ) 
    fprintf ( 1, '\n' );
    fprintf ( 1, 'Error!\n' );
    fprintf ( 1, '  Could not open the output file.\n' );
    error ( 'Error!' );
  end
  for ii = 1 : size(triangle_node,1)
    for jj = 1 : size(triangle_node,2)
      fprintf ( output_unit, '  %12d', triangle_node(ii,jj) );
    end
    fprintf ( output_unit, '\n' );
  end
  fclose ( output_unit );
node_show = 2 ;
triangle_show = 2 ;
triangulation_plot ( 'MG', node_show, triangle_show )






















