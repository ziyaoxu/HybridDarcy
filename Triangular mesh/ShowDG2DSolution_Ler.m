function ShowDG2DSolution_Ler(coordinates,elements3,U)
figure
hold on
quad_lamda1 = [1;0;0];
quad_lamda2 = [0;1;0];
quad_r = 1*quad_lamda1 + (-1)*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
quad_s = (-1)*quad_lamda1 + 1*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
[quad_a,quad_b] = rstoab(quad_r,quad_s) ; 
Np = size(U,1); % degree of freedom in each element
P_order =  (-3+sqrt(1+8*Np))/2 ; % degree of polinomial in each element
basisP = zeros( Np , 3 ) ; % basis function's value in vertex point
for i = 0 : P_order % compute the basis function and it's gradient's value
   for j= 0 : P_order % in quadrature point 
       if (i+j) <= P_order % phi(r,s) = r^i * s^j 
           m = i+(P_order+1)*j+1-(j*(j-1))/2; 
           basisP(m,:) = Simplex2DP(quad_a,quad_b,i,j)' ;  
       end
   end
end
U_3 = zeros(3,size(U,2));
for j=1:size(elements3,2)
    U_3(:,j) = basisP'*U(:,j);
end    
for j=1:size(elements3,2)
    trisurf([1 2 3],coordinates(elements3(:,j),1),coordinates(elements3(:,j),2),U_3(:,j),'facecolor','interp');
end    
hold off
%colormap('jet');
%colorbar('YTickLabel',...
%    {'-1.0','-0.8','-0.6','-0.4','-0.2','0.0','0.2'...
%        '0.4','0.6','0.8','1.0'});           % 
%caxis([-1,1])        % color (saturation)ranges form 0 to 1
%axis image          % 效果与命令axis equal相同(使在每个方向的数据单位都相同)...
% 只是图形区域刚好紧紧包围图象数据。
%axis off            % Remove axes
shading interp % Eliminate the element edge line
end