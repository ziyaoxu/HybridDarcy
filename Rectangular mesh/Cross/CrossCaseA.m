function CrossCaseA % Case (a) : Crossing fracture
% time: 2020.5.28 - 2020.5.28
% author : xuziyao
% Basis function is orthogonal polynomials(any degree) on the rectangular element.
% Dirichlet boundary condition 
% Equation:
% -div(K grad(p)) = f in \Omega
% p = p_D on \Gamma_D, u.n = q_N on \Gamma_N
% =>
% s = -grad p
% u = (km+epsilon*kf*delta*1*nu*nu)s
% q + div(u) = 0 (q=-f)
% =>
% A0*Sx=B1*P+D1, A0*Sy=B2*P+D2, 
% A0*Ux=As11*Sx+As12*Sy
% A0*Uy=As21*Sx+As22*Sy
% A0*Q=C1*Ux+C2*Uy+Nv+Em*P+Ed
% P = (C1*A0inv*F11*A0inv*B1+C1*A0inv*F12*A0inv*B2+C2*A0inv*F21*A0inv*B1+C2*A0inv*F22*A0inv*B2)\
%        (A0Q-Nv-Ed-C1*A0inv*F11*A0inv*D1-C1*A0inv*F12*A0inv*D2-C2*A0inv*F21*A0inv*D1-C2*A0inv*F22*A0inv*D2).
% FINISHED
clc,clear
format long
% set the geometry and triangulation parameters:
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain 
% x direction is divided into N parts ,y direction is divided into M parts:
[Cell_N,Cell_M] = deal(20);
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,~] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );
EtoEmap = Rectangulation_neighbor_rectangles( Cell_N,Cell_M ); % non-periodic boundary condition
Tarea =  hx*hy; Elenth = [hx;hy;hx;hy];
Jacobimat = [2/hx,0;0,2/hy];%Jacobimat=[Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
alpha = 5/hx^3; % penalty on p jump to make the pressure continuous cross interfaces
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb
% ,the # of elements passed by this fracture. 
NumberofFractures = 2; 
ParametersofFracture = zeros(NumberofFractures, 12); 
ParametersofFracture(:,1:4) = [0.5,0.5,0.5,0; 0.5,0.5,0.5,pi/2;]+(rand(NumberofFractures,4))*1e-7; % position of fractures
ParametersofFracture(:,5) = 1e-3; % width of fractures
ParametersofFracture(:,6)= 1e8 ; % tangential permeability_nu of fractures
ParametersofFracture(:,7)= 1e8 ; % normal permeability_nu of fractures
ParametersofFracture(:,8) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,11) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
% plot the triangulation and fractures
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;view([0,90]);hold on
plot(ParametersofFracture(:,[8,10])',ParametersofFracture(:,[9,11])',...
   'r-','LineWidth',2); hold off; drawnow;
% (2) compute basis data
HighOrderDegree = 1 ; % polynomial's order in each element 
HighOrderNp = (HighOrderDegree+1)^2; % degree of freedom in each element for high order polynomial
Gauss_num = 4 ; % the number of quadrature point in [-1,1] 
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
basis1DP = zeros( HighOrderDegree+1 , Gauss_num ) ; % basis function's value in quadrature point
Gradbasis1DP = zeros( HighOrderDegree+1 , Gauss_num ) ; % gradient's value in quadrature point
basis1DPendpoints = zeros( HighOrderDegree+1 , 2 ) ; % basis function's value at -1 and 1
Gradbasis1DPendpoints = zeros( HighOrderDegree+1 , 2 ) ; % gradient's value at -1 and 1
for i = 0 : HighOrderDegree
    basis1DP(i+1,:)     = JacobiP(Gauss_x,0,0,i)';
    Gradbasis1DP(i+1,:) = GradJacobiP(Gauss_x,0,0,i)';
    basis1DPendpoints(i+1,:) = JacobiP([-1;1],0,0,i)';
    Gradbasis1DPendpoints(i+1,:) = GradJacobiP([-1;1],0,0,i)';
end
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(w) = 1
% Tensor product of Lobatto quadrature rule
[quad_lamda1,quad_lamda2]=meshgrid(Gauss_x,Gauss_x); 
quad_lamda1 = quad_lamda1(:);
quad_lamda2 = quad_lamda2(:);
quad_w = Gauss_w*Gauss_w';
quad_w = quad_w(:);
[basisP,GradrbasisP,GradsbasisP] = deal( zeros( HighOrderNp , size(quad_w,1) ) ); % basis function's value in quadrature point
[boundaryP,GradrboundaryP,GradsboundaryP] = deal( zeros(HighOrderNp,Gauss_num,4) );
[x_index,y_index] = meshgrid(1:HighOrderDegree+1,1:HighOrderDegree+1);
x_index = x_index(:); y_index = y_index(:);
for i = 1 : HighOrderNp
    basisP(i,:) = reshape(basis1DP(y_index(i),:)'*basis1DP(x_index(i),:),[],1); 
    GradrbasisP(i,:) = reshape(basis1DP(y_index(i),:)'*Gradbasis1DP(x_index(i),:),[],1); 
    GradsbasisP(i,:) = reshape(Gradbasis1DP(y_index(i),:)'*basis1DP(x_index(i),:),[],1); 
    boundaryP(i,:,1) = reshape(basis1DP(x_index(i),:)* basis1DPendpoints(y_index(i),1),[],1); % lower edge 
    boundaryP(i,:,2) = reshape(basis1DPendpoints(x_index(i),end)* basis1DP(y_index(i),:),[],1); % right edge
    boundaryP(i,:,3) = reshape(basis1DP(x_index(i),end:-1:1)* basis1DPendpoints(y_index(i),end),[],1); % upper edge 
    boundaryP(i,:,4) = reshape(basis1DPendpoints(x_index(i),1)* basis1DP(y_index(i),end:-1:1),[],1); % left edge 
    GradrboundaryP(i,:,1) = reshape(Gradbasis1DP(x_index(i),:)* basis1DPendpoints(y_index(i),1),[],1); 
    GradrboundaryP(i,:,2) = reshape(Gradbasis1DPendpoints(x_index(i),end)* basis1DP(y_index(i),:),[],1); 
    GradrboundaryP(i,:,3) = reshape(Gradbasis1DP(x_index(i),end:-1:1)* basis1DPendpoints(y_index(i),end),[],1); 
    GradrboundaryP(i,:,4) = reshape(Gradbasis1DPendpoints(x_index(i),1)* basis1DP(y_index(i),end:-1:1),[],1); 
    GradsboundaryP(i,:,1) = reshape(basis1DP(x_index(i),:)* Gradbasis1DPendpoints(y_index(i),1),[],1);
    GradsboundaryP(i,:,2) = reshape(basis1DPendpoints(x_index(i),end)* Gradbasis1DP(y_index(i),:),[],1);
    GradsboundaryP(i,:,3) = reshape(basis1DP(x_index(i),end:-1:1)* Gradbasis1DPendpoints(y_index(i),end),[],1);
    GradsboundaryP(i,:,4) = reshape(basis1DPendpoints(x_index(i),1)* Gradbasis1DP(y_index(i),end:-1:1),[],1);
end
% compute fractures basic data
[ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradrbasisPonFracture,GradsbasisPonFracture, ...
    FracturesStartRefCoord,FracturesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements4,ParametersofFracture,Gauss_x,HighOrderDegree);
% (3) define the variables
[HighOrderP,HighOrderUx,HighOrderUy] = deal( zeros( HighOrderNp , size(elements4,2) ) ); % unknowns
[IMat,JMat] = ...
    deal(ones(HighOrderNp*size(elements4,2)*(size(EtoEmap,1)+1),HighOrderNp));
[VA0,VA0inv,VAs11,VAs12,VAs21,VAs22,VB1,VB2,VC1,VC2,VEm] = ...
    deal(zeros(HighOrderNp*size(elements4,2)*(size(EtoEmap,1)+1),HighOrderNp));
[D1,D2,Nv,Ed,A0Q] = deal(zeros(size( elements4,2 )*HighOrderNp,1));
%---------------ELLIPTIC EQUATION LDG SOLVER------------------------------
for jj = 1 : size(elements4,2) % run element from first to end  
xy = [coordinates(elements4(1,jj),1)+hx*quad_lamda1,coordinates(elements4(1,jj),2)+hy*quad_lamda2];
I = (jj-1)*HighOrderNp+(1:HighOrderNp);
II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
[JMat(II,:),IMat(II,:)]   = meshgrid(I,I);
VA0(II,:) = VA0(II,:) + Tarea*basisP*diag(quad_w)*basisP';
VA0inv(II,:) = inv(Tarea*basisP*diag(quad_w)*basisP');
VAs11(II,:)  = VAs11(II,:) + Tarea*basisP*diag(k_m11(xy).*quad_w)*basisP';
VAs12(II,:)  = VAs12(II,:) + Tarea*basisP*diag(k_m12(xy).*quad_w)*basisP';
VAs21(II,:)  = VAs21(II,:) + Tarea*basisP*diag(k_m21(xy).*quad_w)*basisP';
VAs22(II,:)  = VAs22(II,:) + Tarea*basisP*diag(k_m22(xy).*quad_w)*basisP';
VB1(II,:) = VB1(II,:) + Tarea*GradrbasisP*Jacobimat(1,1)*diag(quad_w)*basisP';
VB2(II,:) = VB2(II,:) + Tarea*GradsbasisP*Jacobimat(2,2)*diag(quad_w)*basisP';
VC1(II,:) = VC1(II,:) + Tarea*GradrbasisP*Jacobimat(1,1)*diag(quad_w)*basisP';
VC2(II,:) = VC2(II,:) + Tarea*GradsbasisP*Jacobimat(2,2)*diag(quad_w)*basisP';
A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) = A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
    Tarea*basisP*(q(xy).*quad_w);
% Edge integral
    for ii = 1 : 4  % the ii-th edge of jj-th element
        vertex1 = ii ; 
        vertex2 = mod(ii,4)+1 ;
        edge_vec =  coordinates(elements4(vertex2,jj),:)' - ...
            coordinates(elements4(vertex1,jj),:)' ;
        outer_vec = - [0,-1;1,0]*(edge_vec/norm(edge_vec));
        E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj 
        Boundary_xy = (1-Gauss_x)*coordinates(elements4(vertex1,jj),:) + Gauss_x*coordinates(elements4(vertex2,jj),:);
        if E ~= -1 % the edge ii is not in boundary
            % find the index of the edge in neighborhood element:
            pp =  EtoEmap(:,E) == jj ;   
            J = (E-1)*HighOrderNp+1:E*HighOrderNp ;
            IJ = ((jj-1)*(size(EtoEmap,1)+1)+ii)*HighOrderNp+(1:HighOrderNp);
            [JMat(IJ,:),IMat(IJ,:)] = meshgrid(J,I);
            % central flux p^=(p1+p2)/2, u^=(u1+u2)/2
            VB1(II,:) = VB1(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB1(IJ,:) = VB1(IJ,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VB2(II,:) = VB2(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(IJ,:) = VB2(IJ,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC1(II,:) = VC1(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC1(IJ,:) = VC1(IJ,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC2(II,:) = VC2(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(IJ,:) = VC2(IJ,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            % penalty on jump of p:
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VEm(IJ,:) = VEm(IJ,:) + alpha*Elenth(ii)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
        elseif abs(outer_vec(1))>1/2 % Dirichlet boundary.
            D1((jj-1)*HighOrderNp+(1:HighOrderNp)) = D1((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*(Gauss_w.*p(Boundary_xy));
            D2((jj-1)*HighOrderNp+(1:HighOrderNp)) = D2((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*(Gauss_w.*p(Boundary_xy));
            VC1(II,:) = VC1(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(II,:) = VC2(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) = Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
                alpha*Elenth(ii)*boundaryP(:,:,ii)*(Gauss_w.*p(Boundary_xy));
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
        else % Neumann boundary
            VB1(II,:) = VB1(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(II,:) = VB2(II,:) - Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) = Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii)*boundaryP(:,:,ii)*(Gauss_w.*q_N(Boundary_xy,outer_vec));
        end
    end
end
% proceed the fractures
for k = 1 : NumberofFractures
    FractureNu = [cos(ParametersofFracture(k,4));sin(ParametersofFracture(k,4))];
    FractureSigma = [cos(pi/2+ParametersofFracture(k,4));sin(pi/2+ParametersofFracture(k,4))];
    for ell = 1 : ParametersofFracture(k,12)
    CurrentIndex = sum(ParametersofFracture(1:k-1,12)) + ell;
    jj = FracturesPath( CurrentIndex );
    II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
    VAs11(II,:) = VAs11(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(1)*FractureNu(1);
    VAs12(II,:) = VAs12(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(1)*FractureNu(2);
    VAs21(II,:) = VAs21(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(2)*FractureNu(1);
    VAs22(II,:) = VAs22(II,:)  + ...
        ParametersofFracture(k,5)*ParametersofFracture(k,6)*FracturesLength( CurrentIndex )*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonFracture((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*FractureNu(2)*FractureNu(2);
    end
end
%A0 = sparse(IMat(:),JMat(:),VA0(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
As11 = sparse(IMat(:),JMat(:),VAs11(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
As12 = sparse(IMat(:),JMat(:),VAs12(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
As21 = sparse(IMat(:),JMat(:),VAs21(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
As22 = sparse(IMat(:),JMat(:),VAs22(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
A0inv = sparse(IMat(:),JMat(:),VA0inv(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
B1 = sparse(IMat(:),JMat(:),VB1(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
B2 = sparse(IMat(:),JMat(:),VB2(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
C1 = sparse(IMat(:),JMat(:),VC1(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
C2 = sparse(IMat(:),JMat(:),VC2(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Em = sparse(IMat(:),JMat(:),VEm(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
HighOrderP(:) = (C1*A0inv*As11*A0inv*B1+C1*A0inv*As12*A0inv*B2+C2*A0inv*As21*A0inv*B1+C2*A0inv*As22*A0inv*B2+Em)...
    \(A0Q-Nv-Ed-C1*A0inv*As11*A0inv*D1-C1*A0inv*As12*A0inv*D2-C2*A0inv*As21*A0inv*D1-C2*A0inv*As22*A0inv*D2);
HighOrderUx(:) = (A0inv*As11*A0inv*B1+A0inv*As12*A0inv*B2)*HighOrderP(:)+A0inv*As11*A0inv*D1+A0inv*As12*A0inv*D2;
HighOrderUy(:) = (A0inv*As21*A0inv*B1+A0inv*As22*A0inv*B2)*HighOrderP(:)+A0inv*As21*A0inv*D1+A0inv*As22*A0inv*D2;
HighOrderUx(2:end,FracturesPath)=0;
HighOrderUy(2:end,FracturesPath)=0;
ShowDG2DPressure_Ler(coordinates,elements4,HighOrderP);
%ShowDG2DPressure_Ler(coordinates,elements4,HighOrderUx);
%ShowDG2DPressure_Ler(coordinates,elements4,HighOrderUy);
shading faceted; view([30,30])
%--------------------------------------------------------------------------
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = 1-x;
end
function f = f(xy)
x=xy(:,1);
y=xy(:,2);
f = zeros(size(x));
end
function k_m11 = k_m11(xy)
x=xy(:,1);
y=xy(:,2);
k_m11 = ones(size(x));
end
function k_m12 = k_m12(xy)
x=xy(:,1);
y=xy(:,2);
k_m12 = zeros(size(x));
end
function k_m21 = k_m21(xy)
x=xy(:,1);
y=xy(:,2);
k_m21 = zeros(size(x));
end
function k_m22 = k_m22(xy)
x=xy(:,1);
y=xy(:,2);
k_m22 = ones(size(x));
end
function q_N = q_N(xy,n)
x=xy(:,1);
y=xy(:,2);
q_N = zeros(size(x));
end
function q = q(xy)
q = -f(xy);
end
%--------------------------------------------------------------------------
end