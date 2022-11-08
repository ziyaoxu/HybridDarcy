function [basisP,GradrbasisP,GradsbasisP,boundaryP,GradrboundaryP,GradsboundaryP,vertexP] = ...
    SetLegendreBasisData(P_order,quad_lamda1,quad_lamda2,Gauss_x)
%UNTITLED2 Summary of this function goes here 
% Gauss_x is distributed in [0,1].
% quad_lamda1 , quad_lamda2 are distributed in lower left half of [0,1]*[0,1].
Np = (P_order+1) * (P_order+2) / 2; % degree of freedom in each element
quad_num = size( quad_lamda1 , 1 ); % number of quadrature point on triangle
Gauss_num = size( Gauss_x , 1 );  % number of quadrature point on edge.
basisP = zeros( Np , quad_num ) ; % basis function's value in quadrature point
GradrbasisP = zeros( Np , quad_num ) ; % gradient's value in quadrature point
GradsbasisP = zeros( Np , quad_num ) ; % gradient's value in quadrature point
quad_r = 1*quad_lamda1 + (-1)*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
quad_s = (-1)*quad_lamda1 + 1*quad_lamda2 + (-1)*(1-quad_lamda1-quad_lamda2);
[quad_a,quad_b] = rstoab(quad_r,quad_s) ; 
for ii = 0 : P_order % compute the basis function and it's gradient's value
   for jj = 0 : P_order % in quadrature point 
       if (ii+jj) <= P_order % phi(r,s) = r^i * s^j 
           m = ii+(P_order+1)*jj+1-(jj*(jj-1))/2; 
           basisP(m,:) = Simplex2DP(quad_a,quad_b,ii,jj)' ;  
           [dmodedr, dmodeds] = GradSimplex2DP(quad_a,quad_b,ii,jj);
           GradrbasisP(m,:) =  dmodedr' ;  
           GradsbasisP(m,:) =  dmodeds' ;  
       end
   end
end

quad_rs = zeros( Gauss_num,2,3 ) ;
quad_ab = zeros( Gauss_num,2,3 ) ; 
% coordinate(rs) of Gauss points of edge A1A2 in reference triangular.
quad_rs(:,:,1) = ones(size(Gauss_x))*[-1,1] + Gauss_x*([-1,-1]-[-1,1]) ; 
% coordinate(rs) of Gauss points of edge A2A3 in reference triangular.
quad_rs(:,:,2) = ones(size(Gauss_x))*[-1,-1]+ Gauss_x*([1,-1]-[-1,-1]) ; 
% coordinate(rs) of Gauss points of edge A3A1 in reference triangular.
quad_rs(:,:,3) = ones(size(Gauss_x))*[1,-1] + Gauss_x*([-1,1]-[1,-1])  ; 
for kk = 1 : 3
    [quad_ab(:,1,kk),quad_ab(:,2,kk)] = rstoab(quad_rs(:,1,kk),quad_rs(:,2,kk)) ; 
end
boundaryP = zeros(Np,Gauss_num,3);
GradrboundaryP= zeros(Np,Gauss_num,3);
GradsboundaryP= zeros(Np,Gauss_num,3);
for ii = 0 : P_order 
   for jj= 0 : P_order  
       if (ii+jj) <= P_order % phi(r,s) = r^i * s^j 
           m = ii+(P_order+1)*jj+1-(jj*(jj-1))/2 ; 
           for kk = 1 : 3
           boundaryP(m,:,kk) = Simplex2DP(quad_ab(:,1,kk),quad_ab(:,2,kk),ii,jj)' ;  
           [dmodedr, dmodeds] = GradSimplex2DP(quad_ab(:,1,kk),quad_ab(:,2,kk),ii,jj);
           GradrboundaryP(m,:,kk) = dmodedr' ;  
           GradsboundaryP(m,:,kk) = dmodeds' ;  
           end
       end
   end
end

vertexP = zeros( Np , 3 ) ; % basis function's value in vertex point
quad_vr = [1;-1;-1];
quad_vs = [-1;1;-1];
[quad_va,quad_vb] = rstoab(quad_vr,quad_vs) ; 
for ii = 0 : P_order % compute the basis function and it's gradient's value
   for jj= 0 : P_order % in quadrature point 
       if (ii+jj) <= P_order % phi(r,s) = r^i * s^j 
           m = ii+(P_order+1)*jj+1-(jj*(jj-1))/2; 
           vertexP(m,:) = Simplex2DP(quad_va,quad_vb,ii,jj)' ;  
       end
   end
end

end

