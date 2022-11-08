function h = ShowDG2DSolution_Ler(coordinates,elements4,U)
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
%colormap('parula');
% Display a colorbar beside the axes and use descriptive textstrings as y-tick labels.
colorbar('YTickLabel',...
    {'0.0','0.1','0.2','0.3','0.4','0.5','0.6'...
        '0.7','0.8','0.9','1.0'});           % 
caxis([0,1])        % color (saturation)ranges form 0 to 1
axis image          % Ч��������axis equal��ͬ(ʹ��ÿ���������ݵ�λ����ͬ)...
% ֻ��ͼ������պý�����Χͼ����ݡ�
%axis off            % Remove axes
shading interp % Eliminate the element edge line
end
