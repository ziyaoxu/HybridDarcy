function [ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradrbasisPonFracture,GradsbasisPonFracture, ...
    FracturesStartRefCoord,FracturesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements4,ParametersofFracture,quad_x,P_order)
% time : 2019.7.27 - 2019.7.28
% author : xuziyao
% quad_x is distributed in [0,1].
% Basis function is the orthogonal polynomials(any degree) on the rectangular element.
HighOrderNp = (P_order+1)^2; % degree of freedom in each element for high order polynomial
ParametersofFracture(:,8) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,11) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,12) = 0; % the number of elements which is passed by the k-th fracture
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_NM = size(elements4,2); Cell_N = elements4(4,1)-elements4(1,1)-1; Cell_M = Cell_NM/Cell_N;
EtoEmap = Rectangulation_neighbor_rectangles( Cell_N,Cell_M );
hx = ( max(coordinates(:,1)) - min(coordinates(:,1)) )/ Cell_N ; 
hy = ( max(coordinates(:,2)) - min(coordinates(:,2)) )/ Cell_M ; 
NumberofFractures = size(ParametersofFracture,1);
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb,
% compute the # of elements which are passed through by this fracture. 
StartEndElements = zeros(NumberofFractures,2);
for jj = 1 : size(elements4,2) % compute the start and end elements of each fracture's path
    isStartElement = inpolygon( ParametersofFracture(:,8),ParametersofFracture(:,9),...
        coordinates(elements4([1;2;3;4;1],jj),1),coordinates(elements4([1;2;3;4;1],jj),2) );
    isEndElement   = inpolygon( ParametersofFracture(:,10),ParametersofFracture(:,11),...
        coordinates(elements4([1;2;3;4;1],jj),1),coordinates(elements4([1;2;3;4;1],jj),2) );
    StartEndElements(isStartElement,1) = jj;
    StartEndElements(isEndElement,2) = jj;
end
for k = 1 : NumberofFractures % compute ParametersofFracture(:,12), the number of elements along each fracture's path
    CurrentElement = StartEndElements(k,1); 
    ParametersofFracture(k,12) = ParametersofFracture(k,12) + 1;
    pp = 0;
    while ( CurrentElement ~= StartEndElements(k,2) )
        [~, ~, ii_int] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
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
[x_index,y_index] = meshgrid(1:(P_order+1),1:(P_order+1));
x_index = x_index(:); y_index = y_index(:);
basisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
GradrbasisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
GradsbasisPonFracture = zeros( HighOrderNp*length(FracturesPath) , length(quad_x) + 2 );
FracturesStartRefCoord = zeros(length(FracturesPath),2);
FracturesEndRefCoord = zeros(length(FracturesPath),2);
CurrentIndex = 0;
for k = 1 : NumberofFractures
% compute the basis data for each element which is passed through by the
% k-th fracture:
    CurrentElement = StartEndElements(k,1); CurrentIndex = CurrentIndex + 1;
    FracturesPath(CurrentIndex) = CurrentElement;
    pp = 0; pplamda = ([ParametersofFracture(k,8);ParametersofFracture(k,9)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
    while ( CurrentElement ~= StartEndElements(k,2) )
        [x_int, y_int, ii_int] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
           ,ParametersofFracture(k,[8,10]),ParametersofFracture(k,[9,11])); edge_int = ii_int(:,1);
        qq = setdiff(edge_int,pp); qqlamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
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
        basisP_X = zeros( P_order+1 , length(quad_lamda1) ) ; % basis function's value in quadrature point
        GradbasisP_X = zeros( P_order+1 , length(quad_lamda1) ) ; % gradient's value in quadrature point
        basisP_Y = zeros( P_order+1 , length(quad_lamda2) ) ; % basis function's value in quadrature point
        GradbasisP_Y = zeros( P_order+1 , length(quad_lamda2) ) ; % gradient's value in quadrature point
        for i = 0 : P_order
            basisP_X(i+1,:)     = JacobiP(quad_lamda1*2-1,0,0,i)';
            GradbasisP_X(i+1,:) = GradJacobiP(quad_lamda1*2-1,0,0,i)';
            basisP_Y(i+1,:)     = JacobiP(quad_lamda2*2-1,0,0,i)';
            GradbasisP_Y(i+1,:) = GradJacobiP(quad_lamda2*2-1,0,0,i)';
        end
        for i = 1 : HighOrderNp
            basisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = basisP_X(x_index(i),:).*basisP_Y(y_index(i),:); 
            GradrbasisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = GradbasisP_X(x_index(i),:).*basisP_Y(y_index(i),:); 
            GradsbasisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = basisP_X(x_index(i),:).*GradbasisP_Y(y_index(i),:); 
        end
        NextElement = EtoEmap(qq,CurrentElement);
        pp = find(EtoEmap(:,NextElement) == CurrentElement); pplamda = ([x_int(edge_int==qq);y_int(edge_int==qq)]-coordinates(elements4(1,NextElement),:)')./[hx;hy];
        CurrentElement = NextElement; CurrentIndex = CurrentIndex + 1;
        FracturesPath(CurrentIndex) = CurrentElement;
    end
    qqlamda = ([ParametersofFracture(k,10);ParametersofFracture(k,11)]-coordinates(elements4(1,CurrentElement),:)')./[hx;hy];
    if CurrentElement == StartEndElements(k,1)
    FracturesLength(CurrentIndex) = ...
           sqrt( (ParametersofFracture(k,8)-ParametersofFracture(k,10))^2 + (ParametersofFracture(k,9)-ParametersofFracture(k,11))^2 );    
    else
    [x_int, y_int, ~] = polyxpoly(coordinates(elements4([1;2;3;4;1],CurrentElement),1),coordinates(elements4([1;2;3;4;1],CurrentElement),2)...
       ,ParametersofFracture(k,[8,10]),ParametersofFracture(k,[9,11]));    
    FracturesLength(CurrentIndex) = ...
           sqrt( (x_int-ParametersofFracture(k,10))^2 + (y_int-ParametersofFracture(k,11))^2 );
    end
        % compute basis data (basisPonFracture,GradnubasisPonFracture) in the following:
        quad_lamda1 = pplamda(1) + (qqlamda(1)-pplamda(1))*[quad_x;0;1] ; 
        quad_lamda2 = pplamda(2) + (qqlamda(2)-pplamda(2))*[quad_x;0;1] ;
        FracturesStartRefCoord(CurrentIndex,:) = pplamda';
        FracturesEndRefCoord(CurrentIndex,:) = qqlamda';
        basisP_X = zeros( P_order+1 , length(quad_lamda1) ) ; % basis function's value in quadrature point
        GradbasisP_X = zeros( P_order+1 , length(quad_lamda1) ) ; % gradient's value in quadrature point
        basisP_Y = zeros( P_order+1 , length(quad_lamda2) ) ; % basis function's value in quadrature point
        GradbasisP_Y = zeros( P_order+1 , length(quad_lamda2) ) ; % gradient's value in quadrature point
        for i = 0 : P_order
            basisP_X(i+1,:)     = JacobiP(quad_lamda1*2-1,0,0,i)';
            GradbasisP_X(i+1,:) = GradJacobiP(quad_lamda1*2-1,0,0,i)';
            basisP_Y(i+1,:)     = JacobiP(quad_lamda2*2-1,0,0,i)';
            GradbasisP_Y(i+1,:) = GradJacobiP(quad_lamda2*2-1,0,0,i)';
        end
        for i = 1 : HighOrderNp
            basisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = basisP_X(x_index(i),:).*basisP_Y(y_index(i),:); 
            GradrbasisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = GradbasisP_X(x_index(i),:).*basisP_Y(y_index(i),:); 
            GradsbasisPonFracture((CurrentIndex-1)*HighOrderNp+i,:) = basisP_X(x_index(i),:).*GradbasisP_Y(y_index(i),:); 
        end
end
end