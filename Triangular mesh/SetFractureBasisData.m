function [ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture,...
    FracturesStartRefCoord,FracturesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofFracture,quad_x,P_order)
% time : 2020.6.4 - 2020.6.4
% author : xuziyao
% quad_x is distributed in [0,1].
% Basis function is orthogonal polynomials(any degree) on the triangular element.
ParametersofFracture(:,12) = 0; % the number of elements which is passed by the k-th fracture
EtoEmap = triangulation_neighbor_triangles ( elements3 ); 
NumberofFractures = size(ParametersofFracture,1);
HighOrderNp = (P_order+1) * (P_order+2) / 2; % degree of freedom in each element
% ParametersofFracture: x_c,y_c,length,theta,width,permeability,xa,ya,xb,yb,
% # of elements which are passed through by this fracture. 
ParametersofFracture(:,8) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,11) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
StartEndElements = zeros(NumberofFractures,2);
for jj = 1 : size(elements3,2) % compute the start and end elements of each fracture's path
    isStartElement = inpolygon( ParametersofFracture(:,8),ParametersofFracture(:,9),...
        coordinates(elements3([2;3;1;2],jj),1),coordinates(elements3([2;3;1;2],jj),2) );
    isEndElement   = inpolygon( ParametersofFracture(:,10),ParametersofFracture(:,11),...
        coordinates(elements3([2;3;1;2],jj),1),coordinates(elements3([2;3;1;2],jj),2) );
    StartEndElements(isStartElement,1) = jj;
    StartEndElements(isEndElement,2) = jj;
end
for k = 1 : NumberofFractures % compute ParametersofFracture(:,12), the number of elements along each fracture's path
    CurrentElement = StartEndElements(k,1); 
    ParametersofFracture(k,12) = ParametersofFracture(k,12) + 1;
    pp = 0;
    while ( CurrentElement ~= StartEndElements(k,2) )
        [~, ~, ii_int] = polyxpoly(coordinates(elements3([2;3;1;2],CurrentElement),1),coordinates(elements3([2;3;1;2],CurrentElement),2)...
           ,ParametersofFracture(k,[8,10]),ParametersofFracture(k,[9,11])); edge_int = ii_int(:,1);
        qq = setdiff(edge_int,pp);
        NextElement = EtoEmap(qq,CurrentElement);
        pp = find(EtoEmap(:,NextElement) == CurrentElement) ; 
        CurrentElement = NextElement;
        ParametersofFracture(k,12) = ParametersofFracture(k,12) + 1;
    end
end
FracturesPath = zeros(sum(ParametersofFracture(:,12)),1);
FracturesLength = zeros(sum(ParametersofFracture(:,12)),1);
basisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
GradnubasisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
GradsigmabasisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
FracturesStartRefCoord = zeros(length(FracturesPath),3);
FracturesEndRefCoord = zeros(length(FracturesPath),3);
CurrentIndex = 0;
for k = 1 : NumberofFractures
% compute the basis data for each element which is passed through by the
% k-th fracture:
    CurrentElement = StartEndElements(k,1); CurrentIndex = CurrentIndex + 1;
    FracturesPath(CurrentIndex) = CurrentElement;
    pp = 0; pplamda = [1,1,1;coordinates(elements3(:,CurrentElement),:)'] \ ...
              [1;ParametersofFracture(k,8);ParametersofFracture(k,9)];
    while ( CurrentElement ~= StartEndElements(k,2) )
        [x_int, y_int, ii_int] = polyxpoly(coordinates(elements3([2;3;1;2],CurrentElement),1),coordinates(elements3([2;3;1;2],CurrentElement),2)...
           ,ParametersofFracture(k,[8,10]),ParametersofFracture(k,[9,11])); edge_int = ii_int(:,1);
        qq = setdiff(edge_int,pp); qqlamda = [1,1,1;coordinates(elements3(:,CurrentElement),:)'] \ ...
              [1;x_int(edge_int==qq);y_int(edge_int==qq)];
          if length(edge_int) == 1 
                FracturesLength(CurrentIndex) = ...
                    sqrt( (x_int-ParametersofFracture(k,8))^2 + (y_int-ParametersofFracture(k,9))^2 );
          elseif length(edge_int) == 2 
                FracturesLength(CurrentIndex) = sqrt( diff(x_int)^2 + diff(y_int)^2 );
          else
              error('intersection points is not one or two');
          end
        % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
        quad_lamda1 = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x;0;1] ; 
        quad_lamda2 = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x;0;1] ; 
        FracturesStartRefCoord(CurrentIndex,:) = pplamda';
        FracturesEndRefCoord(CurrentIndex,:) = qqlamda';
        quad_r = 1*quad_lamda1 + (-1)*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
        quad_s = (-1)*quad_lamda1 + 1*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
        [quad_a,quad_b] = rstoab(quad_r,quad_s) ; 
        for ii = 0 : P_order % compute the basis function and it's gradient's value
        for jj = 0 : P_order % in quadrature point 
            if (ii+jj) <= P_order % phi(r,s) = r^i * s^j 
            m = ii+(P_order+1)*jj+1-(jj*(jj-1))/2; 
            basisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = ...
                Simplex2DP(quad_a,quad_b,ii,jj)' ;  
            [dmodedr, dmodeds] = GradSimplex2DP(quad_a,quad_b,ii,jj);
            GradnubasisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = (...
                [dmodedr,dmodeds] * ...
            ( [1,-1,-1;-1,1,-1]*([1,1,1;coordinates(elements3(:,CurrentElement),:)']\[0,0;1,0;0,1]) * ...
              [cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))] ) )';   
            GradsigmabasisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = (...
                [dmodedr,dmodeds] * ...
            ( [1,-1,-1;-1,1,-1]*([1,1,1;coordinates(elements3(:,CurrentElement),:)']\[0,0;1,0;0,1]) * ...
              [-sin(ParametersofFracture(k,4));cos(ParametersofFracture(k,4))] ) )';   
            end
        end
        end
        NextElement = EtoEmap(qq,CurrentElement);
        pp = find(EtoEmap(:,NextElement) == CurrentElement); pplamda = [1,1,1;coordinates(elements3(:,NextElement),:)'] \ ...
              [1;x_int(edge_int==qq);y_int(edge_int==qq)];
        CurrentElement = NextElement; CurrentIndex = CurrentIndex + 1;
        FracturesPath(CurrentIndex) = CurrentElement;
    end
    [x_int, y_int, ~] = polyxpoly(coordinates(elements3([2;3;1;2],CurrentElement),1),coordinates(elements3([2;3;1;2],CurrentElement),2)...
       ,ParametersofFracture(k,[8,10]),ParametersofFracture(k,[9,11])); 
    qqlamda = [1,1,1;coordinates(elements3(:,CurrentElement),:)'] \ ...
          [1;ParametersofFracture(k,10);ParametersofFracture(k,11)];
    FracturesLength(CurrentIndex) = ...
           sqrt( (x_int-ParametersofFracture(k,10))^2 + (y_int-ParametersofFracture(k,11))^2 );
        % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
        quad_lamda1 = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x;0;1] ; 
        quad_lamda2 = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x;0;1] ; 
        FracturesStartRefCoord(CurrentIndex,:) = pplamda';
        FracturesEndRefCoord(CurrentIndex,:) = qqlamda';
        quad_r = 1*quad_lamda1 + (-1)*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
        quad_s = (-1)*quad_lamda1 + 1*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
        [quad_a,quad_b] = rstoab(quad_r,quad_s) ; 
        for ii = 0 : P_order % compute the basis function and it's gradient's value
        for jj = 0 : P_order % in quadrature point 
            if (ii+jj) <= P_order % phi(r,s) = r^i * s^j 
            m = ii+(P_order+1)*jj+1-(jj*(jj-1))/2; 
            basisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = ...
                Simplex2DP(quad_a,quad_b,ii,jj)' ;  
            [dmodedr, dmodeds] = GradSimplex2DP(quad_a,quad_b,ii,jj);
            GradnubasisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = (...
                [dmodedr,dmodeds] * ...
            ( [1,-1,-1;-1,1,-1]*([1,1,1;coordinates(elements3(:,CurrentElement),:)']\[0,0;1,0;0,1]) * ...
              [cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))] ) )';   
            GradsigmabasisPonFracture( (CurrentIndex-1)*HighOrderNp + m , :) = (...
                [dmodedr,dmodeds] * ...
            ( [1,-1,-1;-1,1,-1]*([1,1,1;coordinates(elements3(:,CurrentElement),:)']\[0,0;1,0;0,1]) * ...
              [-sin(ParametersofFracture(k,4));cos(ParametersofFracture(k,4))] ) )';   
            end
        end
        end
end
end