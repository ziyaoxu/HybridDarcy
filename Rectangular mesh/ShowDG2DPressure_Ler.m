function h = ShowDG2DPressure_Ler(coordinates,elements4,U)
h = figure ; 
hold on
Np = size(U,1) ;P_order = sqrt(Np) - 1 ;
[x_index,y_index] = meshgrid(1:P_order+1,1:P_order+1);
x_index = x_index(:); y_index = y_index(:);
quad_lamda1 = [-1;1;1;-1]; quad_lamda2 = [-1;-1;1;1];
basisP = zeros( Np , 4 ) ; % basis function's value in quadrature point
for i = 1 : Np
    basisP(i,:) = (JacobiP(quad_lamda1,0,0,x_index(i)-1).*JacobiP(quad_lamda2,0,0,y_index(i)-1))';
end
U_4 = zeros(4,size(U,2));
for j=1:size(elements4,2)
    U_4(:,j) = basisP'*U(:,j);
end    
for j=1:size(elements4,2)
    trisurf([1 2 3 4],coordinates(elements4(:,j),1),coordinates(elements4(:,j),2),U_4(:,j),'facecolor','interp');
%    trisurf([1 2 3],coordinates(elements4(:,j),1),coordinates(elements4(:,j),2),U_3(:,j),'facecolor','interp');
end    
hold off
colormap('jet');
shading interp % Eliminate the element edge line
end
