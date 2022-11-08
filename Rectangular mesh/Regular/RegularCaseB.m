function RegularCaseB % only draw slices of reference solution and conforming solution.
% time : 2020.5.28 - 2020.5.28
% author : xuziyao
% Basis function is orthogonal polynomials(any degree) on the rectangular element.
% Dirichlet boundary condition 
% Equation:
% s = -grad p
% [a11,a12;a21,a22]u = s
% q + div(u) = 0 (q=-f)
% =>
% A0*Sx=B1*P+D1+Cp11*Ux+Cp12*Uy, A0*Sy=B2*P+D2+Cp21*Ux+Cp22*Uy, 
% Au11*Ux+Au12*Uy=A0*Sx
% Au21*Ux+Au22*Uy=A0*Sy
% A0*Q=C1*Ux+C2*Uy+Nv+Em*P+Ed
% =>Au=[Au11,Au12;Au21,Au22].
% P=([C1,C2]*Au^{-1}*[B1;B2]+Em)^{-1}*(A0Q-Nv-Ed-[C1,C2]*Au^{-1}*[D1;D2])
% where, 
% D1, D2 from Dirichlet boundary condition of p
% Cp11, Cp12, Cp21, Cp22 from the panelty of jump of flux in p^
% Nv from the Nuemann condition of u
% Em from the panelty of jump of pressure in u^
% Ed from the panelty of jump of pressure (dirichlet condition of p) in u^
% FINISHED
clc,clear
format long
% set the geometry and triangulation parameters :
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts:
[Cell_N,Cell_M] = deal(16);
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,~] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );
EtoEmap = Rectangulation_neighbor_rectangles( Cell_N,Cell_M ); % non-periodic boundary condition
Tarea =  hx*hy; Elenth = [hx;hy;hx;hy];
Jacobimat = [2/hx,0;0,2/hy];%Jacobimat=[Dr/Dx,Dr/Dy ; Ds/Dx,Ds/Dy] 
alpha = 1/hx^2; % penalty on p jump to make the pressure continuous cross interfaces
beta = 5/hx^2; % penalty on u jump to make the flux continuous cross interfaces
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb
% ,the # of elements passed by this fracture. 
NumberofBarriers = 6; NumberofFractures = 0;
ParametersofBarrier = zeros(NumberofBarriers, 12); 
ParametersofBarrier(:,1:4) = [ 0.5-1e-9, 0.5-2e-9, 1-4e-8, 0;...%... % position of barriers
                                0.5-1e-9, 0.5-2e-9, 1-4e-8, pi/2;
                                0.75-1e-9, 0.75-2e-9, 0.5-4e-8, 0;
                                0.75-1e-9, 0.75-2e-9, 0.5-4e-8, pi/2;
                                0.625-1e-9, 0.625-2e-9, 0.25-4e-8, 0;
                                0.625-1e-9, 0.625-2e-9, 0.25-4e-8, pi/2];
ParametersofBarrier(:,5) = 1e-4; % width of fractures
ParametersofBarrier(:,6)= 1e-4 ; % tangential permeability_nu of fractures
ParametersofBarrier(:,7)= 1e-4 ; % normal permeability_nu of fractures
ParametersofBarrier(:,8) = ParametersofBarrier(:,1) - 1/2*ParametersofBarrier(:,3).*cos(ParametersofBarrier(:,4));
ParametersofBarrier(:,9) = ParametersofBarrier(:,2) - 1/2*ParametersofBarrier(:,3).*sin(ParametersofBarrier(:,4));
ParametersofBarrier(:,10) = ParametersofBarrier(:,1) + 1/2*ParametersofBarrier(:,3).*cos(ParametersofBarrier(:,4));
ParametersofBarrier(:,11) = ParametersofBarrier(:,2) + 1/2*ParametersofBarrier(:,3).*sin(ParametersofBarrier(:,4));
% plot the triangulation and fractures
figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;view([0,90]);hold on
plot(ParametersofBarrier(:,[8,10])',ParametersofBarrier(:,[9,11])',...
   'b-','LineWidth',2); hold off; drawnow;
% (2) compute basis data
HighOrderDegree = 1 ; % polynomial's order in each element 
HighOrderNp = (HighOrderDegree+1)^2; % degree of freedom in each element for high order polynomial
Gauss_num = 6 ; % the number of quadrature point in [-1,1] 
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
NumofCorners = 9;
BasisPonCorners = zeros(HighOrderNp,NumofCorners^2); % basis values at NumofCorners^2 points of a rectangle
for i = 1 : HighOrderNp
    for j = 1 : NumofCorners
        for k = 1 : NumofCorners
            BasisPonCorners(i,(j-1)*NumofCorners+k) = JacobiP(-1+2/(NumofCorners-1)*(k-1),0,0,x_index(i)-1) ...
                *JacobiP(-1+2/(NumofCorners-1)*(j-1),0,0,y_index(i)-1);
        end
    end
end
% compute fractures and Barriers basic data
[ParametersofBarrier,BarriersPath,BarriersLength,basisPonBarrier,GradrbasisPonBarrier,GradsbasisPonBarrier, ...
    BarriersStartRefCoord,BarriersEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements4,ParametersofBarrier,Gauss_x,HighOrderDegree);
FluxPenaltyEdges = ...
    JumpPenaltyEdges(EtoEmap,BarriersPath,BarriersStartRefCoord,BarriersEndRefCoord);
% (3) define the variables
[HighOrderP,HighOrderUx,HighOrderUy] = deal( zeros( HighOrderNp , size(elements4,2) ) ); % unknowns
[IMat,JMat] = ...
    deal(ones(HighOrderNp*size(elements4,2)*(size(EtoEmap,1)+1),HighOrderNp));
[VAu11,VAu12,VAu21,VAu22,VCp11,VCp12,VCp21,VCp22,VB1,VB2,VC1,VC2,VEm] = ...
    deal(zeros(HighOrderNp*size(elements4,2)*(size(EtoEmap,1)+1),HighOrderNp));
[D1,D2,Nv,Ed,A0Q] = deal(zeros(size( elements4,2 )*HighOrderNp,1));
%---------------ELLIPTIC EQUATION LDG SOLVER------------------------------
for jj = 1 : size(elements4,2) % run element from first to end  
xy = [coordinates(elements4(1,jj),1)+hx*quad_lamda1,coordinates(elements4(1,jj),2)+hy*quad_lamda2];
I = (jj-1)*HighOrderNp+(1:HighOrderNp);
II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
[JMat(II,:),IMat(II,:)]   = meshgrid(I,I);
%VA0(II,:) = VA0(II,:) + Tarea*basisP*diag(quad_w)*basisP';
VAu11(II,:)  = VAu11(II,:) + Tarea*basisP*diag(a11(xy).*quad_w)*basisP';
VAu12(II,:)  = VAu12(II,:) + Tarea*basisP*diag(a12(xy).*quad_w)*basisP';
VAu21(II,:)  = VAu21(II,:) + Tarea*basisP*diag(a21(xy).*quad_w)*basisP';
VAu22(II,:)  = VAu22(II,:) + Tarea*basisP*diag(a22(xy).*quad_w)*basisP';
%VAinv(II,:)  = VAinv(II,:)  + inv(Tarea*basisP*diag(a(xy).*quad_w)*basisP');
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
            if FluxPenaltyEdges(ii,jj) == 1 % penalty on jump of u:
            VCp11(II,:) = VCp11(II,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(1));
            VCp11(IJ,:) = VCp11(IJ,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(1));
            VCp12(II,:) = VCp12(II,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(2));
            VCp12(IJ,:) = VCp12(IJ,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(2));
            VCp21(II,:) = VCp21(II,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(1));
            VCp21(IJ,:) = VCp21(IJ,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(1));
            VCp22(II,:) = VCp22(II,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)'*(outer_vec(2));
            VCp22(IJ,:) = VCp22(IJ,:) - beta*Elenth(ii)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)'*(-outer_vec(2));
            else % penalty on jump of p:
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VEm(IJ,:) = VEm(IJ,:) + alpha*Elenth(ii)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            end
        elseif (outer_vec(1))>1/2 % Dirichlet boundary.
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
% process the Barriers
for k = 1 : NumberofBarriers
    BarrierNu = [cos(ParametersofBarrier(k,4));sin(ParametersofBarrier(k,4))];
    BarrierSigma = [cos(pi/2+ParametersofBarrier(k,4));sin(pi/2+ParametersofBarrier(k,4))];
    for ell = 1 : ParametersofBarrier(k,12)
    CurrentIndex = sum(ParametersofBarrier(1:k-1,12)) + ell;
    jj = BarriersPath( CurrentIndex );
    II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
    VAu11(II,:) = VAu11(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(1)*BarrierSigma(1);
    VAu12(II,:) = VAu12(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(1)*BarrierSigma(2);
    VAu21(II,:) = VAu21(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(2)*BarrierSigma(1);
    VAu22(II,:) = VAu22(II,:)  + ...
        ParametersofBarrier(k,5)/ParametersofBarrier(k,7)*BarriersLength( CurrentIndex )*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)*diag(Gauss_w)*...
        basisPonBarrier((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),1:end-2)'*BarrierSigma(2)*BarrierSigma(2);
    end
end
%A0 = sparse(IMat(:),JMat(:),VA0(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Au11 = sparse(IMat(:),JMat(:),VAu11(:)-VCp11(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Au12 = sparse(IMat(:),JMat(:),VAu12(:)-VCp12(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Au21 = sparse(IMat(:),JMat(:),VAu21(:)-VCp21(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Au22 = sparse(IMat(:),JMat(:),VAu22(:)-VCp22(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
%Ainv = sparse(IMat(:),JMat(:),VAinv(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
B1 = sparse(IMat(:),JMat(:),VB1(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
B2 = sparse(IMat(:),JMat(:),VB2(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
C1 = sparse(IMat(:),JMat(:),VC1(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
C2 = sparse(IMat(:),JMat(:),VC2(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
Em = sparse(IMat(:),JMat(:),VEm(:),size(elements4,2)*HighOrderNp,size(elements4,2)*HighOrderNp);
HighOrderP(:) = ([C1,C2]*([Au11,Au12;Au21,Au22]\[B1;B2])+Em)\(A0Q-Nv-Ed-[C1,C2]*([Au11,Au12;Au21,Au22]\[D1;D2]));
HighOrderU = [Au11,Au12;Au21,Au22]\[B1*HighOrderP(:)+D1;B2*HighOrderP(:)+D2];
HighOrderUx(:) = HighOrderU( 1 : HighOrderNp*size(elements4,2) );
HighOrderUy(:) = HighOrderU( HighOrderNp*size(elements4,2)+1 : end );
%HighOrderP(2:end,BarriersPath)=0;
%HighOrderP(2:end,sum(FluxPenaltyEdges,1)>=size(FluxPenaltyEdges,1)-1)=0;

% draw contour
ShowDG2DPressure_Ler(coordinates,elements4,HighOrderP);
colormap('jet');caxis([1.0,3.6]);cb=colorbar;cb.Position = cb.Position.*[1.07,1,1,1];
axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);
hold on; plot3(ParametersofBarrier(:,[8,10])',ParametersofBarrier(:,[9,11])',100*ones(size(ParametersofBarrier(:,[9,11])')),...
   'w-','LineWidth',0.05); hold off; 
drawnow

% profiles on slices
NumberofSlices = 1; 
ParametersofSlice = zeros(NumberofSlices, 12); 
ParametersofSlice(:,8:11)=([0.0,0.1,0.9,1.0]-0.5)*(1-exp(-7))+0.5;
ParametersofSlice(:,1) =  (ParametersofSlice(:,8)+ParametersofSlice(:,10))/2;% position of fractures
ParametersofSlice(:,2) =  (ParametersofSlice(:,9)+ParametersofSlice(:,11))/2;% position of fractures
ParametersofSlice(:,3) =  sqrt((ParametersofSlice(:,8)-ParametersofSlice(:,10)).^2+(ParametersofSlice(:,9)-ParametersofSlice(:,11)).^2);% position of fractures
ParametersofSlice(:,4) =  atan2(ParametersofSlice(:,11)-ParametersofSlice(:,9),ParametersofSlice(:,10)-ParametersofSlice(:,8));% position of fractures
[ParametersofSlice,SlicesPath,SlicesLength,basisPonSlice,GradrbasisPonSlice,GradsbasisPonSlice, ...
    SlicesStartRefCoord,SlicesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements4,ParametersofSlice,Gauss_x,HighOrderDegree);
for k = 1 : NumberofSlices
figure; hold on;
    for ell = 1 : ParametersofSlice(k,12)
        CurrentIndex = sum(ParametersofSlice(1:k-1,12)) + ell;
        jj = SlicesPath( CurrentIndex );
        SlicesStartPhysicCoord = coordinates(elements4(1,jj),:)+[hx,hy].*SlicesStartRefCoord(CurrentIndex,:);
        SlicesEndPhysicCoord = coordinates(elements4(1,jj),:)+[hx,hy].*SlicesEndRefCoord(CurrentIndex,:);
        Slice_XY = (1-[0;Gauss_x;1])*SlicesStartPhysicCoord+[0;Gauss_x;1]*SlicesEndPhysicCoord;
        Slice_Re = sqrt(sum((Slice_XY-ones(size([0;Gauss_x;1]))*ParametersofSlice(k,8:9)).^2,2));
        Slice_Ue = basisPonSlice((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),[end-1,1:end-2,end])'*HighOrderP(:,jj);
        handle1=plot(Slice_Re,Slice_Ue,'b','LineWidth',1.5);  
    end
    hold off;
end
% format the plot
xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('pressure', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 10);
set(gcf, 'PaperPositionMode', 'auto')
fig_pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)+1])
hold off;

%--------------------------------------------------------------------------
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = ones(size(x)); 
end
function f = f(xy)
x=xy(:,1);
y=xy(:,2);
f = zeros(size(x));
end
function k_m = k_m(xy)
x=xy(:,1);
y=xy(:,2);
k_m = ones(size(x));
end
function q_N = q_N(xy,n)
    x=xy(:,1);
    y=xy(:,2);
    q_N = zeros(size(x));
    if n(1)<-1/2
    q_N(:)=-1;
    end
end
function a11 = a11(xy)
a11 = 1./k_m(xy);
end
function a12 = a12(xy)
a12 = zeros(size(xy,1),1);
end
function a21 = a21(xy)
a21 = zeros(size(xy,1),1);
end
function a22 = a22(xy)
a22 = 1./k_m(xy);
end
function q = q(xy)
q = -f(xy);
end
%--------------------------------------------------------------------------
end