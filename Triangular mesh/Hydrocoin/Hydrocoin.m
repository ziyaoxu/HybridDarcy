function Hydrocoin % Hydrocoin: fine non-conforming meshes with piecewise linear element
% time: 2020.6.4 - 2020.6.5
% author : xuziyao
% Basis function is orthogonal polynomials(any degree) on the triangular element.
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
% set the geometry and triangulation manners:
%--------------------------------------------------------------------------
% mesh 1:(Distmesh)--------------------------------------------------------
%pv=[0,150; 400,100; 800,150; 1200,100; 1600,150; 1600,-1000; 0,-1000]; h0=5; geps=.0001*h0;
%[coordinates,elements3]=distmesh2d(@dpoly,@huniform,10,[0,-1000;1600,160],pv,pv);elements3 = elements3';
% mesh2:(existing mesh)----------------------------------------------------
%pv=[0,150; 400,100; 800,150; 1200,100; 1600,150; 1600,-1000; 0,-1000]; h0=5; geps=.0001*h0;
%elements3_1=0; elements3=0; coordinates=0;
%load('boxdfm_hydrocoin.mat')
%load('ccdfm_hydrocoin.mat')
%elements3=delaunayn(coordinates); % List of triangles
%pmid=(coordinates(elements3(:,1),:)+coordinates(elements3(:,2),:)+coordinates(elements3(:,3),:))/3; % Compute centroids
%elements3=elements3(feval(@dpoly,pmid,pv)<-geps,:); % Keep interior triangles
%elements3 = elements3';
% mesh 3:(Matlab mesh, non-conforming)-----------------------------------------------------
geodesc = [2 7 0 400 800 1200 1600 1600 0 150 100 150 100 150 -1000 -1000]';
h = 60;
geometry = decsg(geodesc); % Create the geometry
model = createpde;
geometryFromEdges(model,geometry);
generateMesh(model,'Hmax',h,'GeometricOrder','linear');
coordinates = model.Mesh.Nodes'; elements3 =  model.Mesh.Elements;
%pdeplotcoordinates',elements3,'EdgeColor','black','ElementLabels','off')
bcol=[1,1,1];figure;simpplot(coordinates,elements3',[],bcol);
%--------------------------------------------------------------------------
% mesh 4:(Matlab mesh,conforming, domain decomposition)-----------------------------------------------------
%geodesc = zeros(12,4); % 4 subdomains
%geodesc(1:(2+5*2),1) = [2;5;0;400;14000/13;1000;0;150;100;-7500/13;-1000;-1000];
%geodesc(1:(2+4*2),2) = [2;4;400;800;1200;14000/13;100;150;100;-7500/13];
%geodesc(1:(2+5*2),3) = [2;5;1200;1600;1600;1500;14000/13;100;150;-1000;-1000;-7500/13];
%geodesc(1:(2+3*2),4) = [2;3;14000/13;1500;1000;-7500/13;-1000;-1000];
%ns = char('A','B','C','D'); ns = ns';
%sf = 'A+B+C+D';
%h = 1000;
%geometry = decsg(geodesc,sf,ns); % Create the geometry
%model = createpde;
%geometryFromEdges(model,geometry);
%generateMesh(model,'Hmax',h,'GeometricOrder','linear');
%coordinates = model.Mesh.Nodes'; elements3 =  model.Mesh.Elements;
%pdeplot(coordinates',elements3,'EdgeColor','black','ElementLabels','off');
%bcol=[1,1,1];figure;simpplot(coordinates,elements3',[],bcol);hold on;
%--------------------------------------------------------------------------
EtoEmap = triangulation_neighbor_triangles ( elements3 ); 
[ Tarea , Elenth , Jacobimat ] = ComputeGeometryData( elements3 ,coordinates );
alpha = 1/min(Tarea)^1; % penalty on p jump to make the pressure continuous cross interfaces
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_t,permeability_n,xa,ya,xb,yb,# of elements passed by this fracture. 
NumberofFractures = 2; 
ParametersofFracture = zeros(NumberofFractures, 12); 
ParametersofFracture(:,1:4) = [ 950+1e-8, -450, 1100*1.41421356, -pi/4;... % position of fractures
                                1100+1e-8, -450, 500*2.23606797, atan(11/2)  ];
ParametersofFracture(:,5) = [5*sqrt(2);33*sqrt(5)/5]; % width of fractures
ParametersofFracture(:,6)= 1e-6 ; % permeability of fractures
ParametersofFracture(:,7)= 1e-6 ; % permeability of fractures
% plot the triangulation and fractures
ParametersofFracture(:,8) = ParametersofFracture(:,1) - 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,9) = ParametersofFracture(:,2) - 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
ParametersofFracture(:,10) = ParametersofFracture(:,1) + 1/2*ParametersofFracture(:,3).*cos(ParametersofFracture(:,4));
ParametersofFracture(:,11) = ParametersofFracture(:,2) + 1/2*ParametersofFracture(:,3).*sin(ParametersofFracture(:,4));
%simpplot(coordinates,elements3');hold on
hold on; plot(ParametersofFracture(:,[8,10])',ParametersofFracture(:,[9,11])',...
   'k-','LineWidth',1.5); hold off; drawnow;
% (2) compute basis data
HighOrderDegree = 1 ; % polynomial's order in each element 
HighOrderNp = (HighOrderDegree+1) * (HighOrderDegree+2) / 2; % degree of freedom in each element for high order polynomial
Gauss_num = 6 ; % the number of Gauss quadrature point in [-1,1]
Lobatto_num = 6 ; % the number of Lobatto quadrature point in [-1,1]
[Gauss_x,Gauss_w] = JacobiGQ(0,0,Gauss_num-1); 
Gauss_x = 0 + ( Gauss_x - (-1) ) / 2 ; % Gauss_x now is in [0,1].
Gauss_w = Gauss_w / 2 ;% sum(Gauss_w) = 1
[ quad_lamda1 , quad_lamda2 , quad_w ] = quad_rule1 ( 21 ) ; % positive quadrature rule
[basisP,GradrbasisP,GradsbasisP,boundaryP,GradrboundaryP,GradsboundaryP] = ...
    SetLegendreBasisData(HighOrderDegree,quad_lamda1,quad_lamda2,Gauss_x) ;
% compute fractures and Barriers basic data
[ParametersofFracture,FracturesPath,FracturesLength,basisPonFracture,GradnubasisPonFracture,GradsigmabasisPonFracture,...
    FracturesStartRefCoord,FracturesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofFracture,Gauss_x,HighOrderDegree);
% solver
[HighOrderP,HighOrderUx,HighOrderUy] = PoissonSolver;
ShowDG2DSolution_Ler(coordinates,elements3,HighOrderP)
colormap('jet'); caxis([100,150]);ylabel('z') ;cb=colorbar;cb.Position = cb.Position.*[1.095,1,1,1];
hold on; plot3(ParametersofFracture(:,[8,10])',ParametersofFracture(:,[9,11])',(200)*ones(size(ParametersofFracture(:,[9,11])')),...
   'w-','LineWidth',0.05); hold off; axis equal; axis off;ax=axis;axis(ax*1.001); view([0,90]);

% profiles on slices
NumberofSlices = 1; 
ParametersofSlice = zeros(NumberofSlices, 12); 
ParametersofSlice(:,8:11)=[0+exp(-10),-200,1600-exp(-10),-200];
ParametersofSlice(:,1) =  (ParametersofSlice(:,8)+ParametersofSlice(:,10))/2;% position of fractures
ParametersofSlice(:,2) =  (ParametersofSlice(:,9)+ParametersofSlice(:,11))/2;% position of fractures
ParametersofSlice(:,3) =  sqrt((ParametersofSlice(:,8)-ParametersofSlice(:,10)).^2+(ParametersofSlice(:,9)-ParametersofSlice(:,11)).^2);% position of fractures
ParametersofSlice(:,4) =  atan2(ParametersofSlice(:,11)-ParametersofSlice(:,9),ParametersofSlice(:,10)-ParametersofSlice(:,8));% position of fractures
[ParametersofSlice,SlicesPath,SlicesLength,basisPonSlice,GradrbasisPonSlice,GradsbasisPonSlice, ...
    SlicesStartRefCoord,SlicesEndRefCoord] = ...
    SetFractureBasisData(coordinates,elements3,ParametersofSlice,Gauss_x,HighOrderDegree);
for k = 1 : NumberofSlices
    figure; hold on;
    for ell = 1 : ParametersofSlice(k,12)
        CurrentIndex = sum(ParametersofSlice(1:k-1,12)) + ell;
        jj = SlicesPath( CurrentIndex );
        SlicesStartPhysicCoord = SlicesStartRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:); 
        SlicesEndPhysicCoord = SlicesEndRefCoord(CurrentIndex,:) * coordinates(elements3(:,jj),:); 
        Slice_XY = (1-[0;Gauss_x;1])*SlicesStartPhysicCoord+[0;Gauss_x;1]*SlicesEndPhysicCoord;
        Slice_Re = sqrt(sum((Slice_XY-ones(size([0;Gauss_x;1]))*ParametersofSlice(k,8:9)).^2,2));
        Slice_Ue = basisPonSlice((CurrentIndex-1)*HighOrderNp+(1:HighOrderNp),[end-1,1:end-2,end])'*HighOrderP(:,jj);
        handle1=plot(Slice_Re,Slice_Ue,'b','LineWidth',1.5);
    end    
    hold off;
% format the plot
xlabel('arc length', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
ylabel('solution', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
set(gca, 'FontSize', 10);
set(gcf, 'PaperPositionMode', 'auto')
fig_pos = get(gcf, 'PaperPosition');
set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)+1])
hold off;    
end

function [HighOrderP,HighOrderUx,HighOrderUy] = PoissonSolver
% (3) define the variables
[HighOrderP,HighOrderUx,HighOrderUy] = deal( zeros( HighOrderNp , size(elements3,2) ) ); % unknowns
[IMat,JMat] = ...
    deal(ones(HighOrderNp*size(elements3,2)*(size(EtoEmap,1)+1),HighOrderNp));
[VA0,VA0inv,VAs11,VAs12,VAs21,VAs22,VB1,VB2,VC1,VC2,VEm] = ...
    deal(zeros(HighOrderNp*size(elements3,2)*(size(EtoEmap,1)+1),HighOrderNp));
[D1,D2,Nv,Ed,A0Q] = deal(zeros(size( elements3,2 )*HighOrderNp,1));
%---------------ELLIPTIC EQUATION LDG SOLVER------------------------------
for jj = 1 : size(elements3,2) % run element from first to end  
xy = quad_lamda1 * coordinates(elements3(1,jj),:)... 
     + quad_lamda2 * coordinates(elements3(2,jj),:)... 
     +(1-quad_lamda1-quad_lamda2)*coordinates(elements3(3,jj),:) ; 
I = (jj-1)*HighOrderNp+(1:HighOrderNp);
II = (jj-1)*(size(EtoEmap,1)+1)*HighOrderNp+(1:HighOrderNp);
[JMat(II,:),IMat(II,:)]   = meshgrid(I,I);
VA0(II,:) = VA0(II,:) + Tarea(jj)*basisP*diag(quad_w)*basisP';
VA0inv(II,:) = inv(Tarea(jj)*basisP*diag(quad_w)*basisP');
VAs11(II,:)  = VAs11(II,:) + Tarea(jj)*basisP*diag(k_m11(xy).*quad_w)*basisP';
VAs12(II,:)  = VAs12(II,:) + Tarea(jj)*basisP*diag(k_m12(xy).*quad_w)*basisP';
VAs21(II,:)  = VAs21(II,:) + Tarea(jj)*basisP*diag(k_m21(xy).*quad_w)*basisP';
VAs22(II,:)  = VAs22(II,:) + Tarea(jj)*basisP*diag(k_m22(xy).*quad_w)*basisP';
VB1(II,:) = VB1(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))*diag(quad_w)*basisP';
VB2(II,:) = VB2(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))*diag(quad_w)*basisP';
VC1(II,:) = VC1(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,1,jj)+GradsbasisP*Jacobimat(2,1,jj))*diag(quad_w)*basisP';
VC2(II,:) = VC2(II,:) + Tarea(jj)*(GradrbasisP*Jacobimat(1,2,jj)+GradsbasisP*Jacobimat(2,2,jj))*diag(quad_w)*basisP';
A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) = A0Q((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
    Tarea(jj)*basisP*(q(xy).*quad_w);
% Edge integral
    for ii = 1 : 3  % the ii-th edge of jj-th element
        vertex1 = mod(ii,3)+1 ; 
        vertex2 = mod(ii+1,3)+1 ;
        edge_vec =  coordinates(elements3(vertex2,jj),:)' - ...
            coordinates(elements3(vertex1,jj),:)' ;
        outer_vec = - [0,-1;1,0]*(edge_vec/norm(edge_vec));
        E = EtoEmap(ii,jj) ; % the index of neighborhood of element jj 
        if E ~= -1 % the edge ii is not in boundary
            % find the index of the edge in neighborhood element:
            pp =  EtoEmap(:,E) == jj ;   
            J = (E-1)*HighOrderNp+1:E*HighOrderNp ;
            IJ = ((jj-1)*(size(EtoEmap,1)+1)+ii)*HighOrderNp+(1:HighOrderNp);
            [JMat(IJ,:),IMat(IJ,:)] = meshgrid(J,I);
            % central flux p^=(p1+p2)/2, u^=(u1+u2)/2
            VB1(II,:) = VB1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB1(IJ,:) = VB1(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VB2(II,:) = VB2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(IJ,:) = VB2(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC1(II,:) = VC1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC1(IJ,:) = VC1(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            VC2(II,:) = VC2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(IJ,:) = VC2(IJ,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*1/2*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
            % penalty on jump of p:
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VEm(IJ,:) = VEm(IJ,:) + alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,end:-1:1,pp)';
        elseif outer_vec(2)>1/3 % Dirichlet boundary.
            Boundary_xy = (1-Gauss_x)*coordinates(elements3(vertex1,jj),:) + Gauss_x*coordinates(elements3(vertex2,jj),:);
            D1((jj-1)*HighOrderNp+(1:HighOrderNp)) = D1((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*(Gauss_w.*p(Boundary_xy));
            D2((jj-1)*HighOrderNp+(1:HighOrderNp)) = D2((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*(Gauss_w.*p(Boundary_xy));
            VC1(II,:) = VC1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VC2(II,:) = VC2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) = Ed((jj-1)*HighOrderNp+(1:HighOrderNp)) + ...
                alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*(Gauss_w.*p(Boundary_xy));
            VEm(II,:) = VEm(II,:) - alpha*Elenth(ii,jj)*boundaryP(:,:,ii)*diag(Gauss_w)*boundaryP(:,:,ii)';
        else % Neumann boundary
            Boundary_xy = (1-Gauss_x)*coordinates(elements3(vertex1,jj),:) + Gauss_x*coordinates(elements3(vertex2,jj),:);
            VB1(II,:) = VB1(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(1)*diag(Gauss_w)*boundaryP(:,:,ii)';
            VB2(II,:) = VB2(II,:) - Elenth(ii,jj)*boundaryP(:,:,ii)*outer_vec(2)*diag(Gauss_w)*boundaryP(:,:,ii)';
            Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) = Nv((jj-1)*HighOrderNp+(1:HighOrderNp)) - ...
                Elenth(ii,jj)*boundaryP(:,:,ii)*(Gauss_w.*q_N(Boundary_xy,outer_vec));
        end
    end
end
% process the fractures
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
%A0 = sparse(IMat(:),JMat(:),VA0(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As11 = sparse(IMat(:),JMat(:),VAs11(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As12 = sparse(IMat(:),JMat(:),VAs12(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As21 = sparse(IMat(:),JMat(:),VAs21(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
As22 = sparse(IMat(:),JMat(:),VAs22(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
A0inv = sparse(IMat(:),JMat(:),VA0inv(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
B1 = sparse(IMat(:),JMat(:),VB1(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
B2 = sparse(IMat(:),JMat(:),VB2(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
C1 = sparse(IMat(:),JMat(:),VC1(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
C2 = sparse(IMat(:),JMat(:),VC2(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
Em = sparse(IMat(:),JMat(:),VEm(:),size(elements3,2)*HighOrderNp,size(elements3,2)*HighOrderNp);
HighOrderP(:) = (C1*A0inv*As11*A0inv*B1+C1*A0inv*As12*A0inv*B2+C2*A0inv*As21*A0inv*B1+C2*A0inv*As22*A0inv*B2+Em)...
    \(A0Q-Nv-Ed-C1*A0inv*As11*A0inv*D1-C1*A0inv*As12*A0inv*D2-C2*A0inv*As21*A0inv*D1-C2*A0inv*As22*A0inv*D2);
HighOrderUx(:) = (A0inv*As11*A0inv*B1+A0inv*As12*A0inv*B2)*HighOrderP(:)+A0inv*As11*A0inv*D1+A0inv*As12*A0inv*D2;
HighOrderUy(:) = (A0inv*As21*A0inv*B1+A0inv*As22*A0inv*B2)*HighOrderP(:)+A0inv*As21*A0inv*D1+A0inv*As22*A0inv*D2;
HighOrderUx(2:end,FracturesPath)=0;
HighOrderUy(2:end,FracturesPath)=0;
end
%--------------------------------------------------------------------------
function p = p(xy)
x = xy(:,1);
y = xy(:,2);
p = y; 
end
function f = f(xy)
x=xy(:,1);
y=xy(:,2);
f = zeros(size(x));
end
function k_m11 = k_m11(xy)
x=xy(:,1);
y=xy(:,2);
k_m11 = (1e-8)*ones(size(x)); 
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
k_m22 = (1e-8)*ones(size(x)); 
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