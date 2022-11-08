clc,clear
% set the geometry and triangulation parameters:
xmin = 0 ; xmax = 1 ; % domain 
ymin = 0 ; ymax = 1 ; % domain
% x direction is divided into N parts ,y direction is divided into M parts:
Cell_N = 1 ; Cell_M = 1 ; 
hx = ( xmax - xmin ) / Cell_N ; hy = ( ymax - ymin ) / Cell_M ;
[coordinates,elements4,EtoEmap] = RectangularMesh( xmin,xmax,ymin,ymax,Cell_N,Cell_M );% periodic boundary condition
% plot the triangulation and fractures
% ParametersofFracture: x_c,y_c,length,theta,width,permeability_nu,permeability_sigma,xa,ya,xb,yb
% ,the # of elements passed by this fracture. 
FractureData = [...
 0.0500 0.4160 0.2200 0.0624;
 0.0500 0.2750 0.2500 0.1350;
 0.1500 0.6300 0.4500 0.0900;
 0.7000 0.2350 0.8500 0.1675;
 0.6000 0.3800 0.8500 0.2675;
 0.3500 0.9714 0.8000 0.7143;
 0.7500 0.9574 0.9500 0.8155;
 0.1500 0.8363 0.4000 0.9727];
BarrierData  = [...
 0.1500 0.9167 0.4000 0.5000;
 0.6500 0.8333 0.8500 0.1667];
NumberofFractures = size(FractureData,1); 
ParametersofFracture = zeros(NumberofFractures, 12); 
ParametersofFracture(:,8:11) = FractureData;
ParametersofFracture(:,1) =  (FractureData(:,1)+FractureData(:,3))/2;% position of fractures
ParametersofFracture(:,2) =  (FractureData(:,2)+FractureData(:,4))/2;% position of fractures
ParametersofFracture(:,3) =  sqrt((FractureData(:,1)-FractureData(:,3)).^2+(FractureData(:,2)-FractureData(:,4)).^2);% position of fractures
ParametersofFracture(:,4) =  atan2(FractureData(:,4)-FractureData(:,2),FractureData(:,3)-FractureData(:,1));% position of fractures

NumberofBarriers = size(BarrierData,1); 
ParametersofBarrier = zeros(NumberofBarriers, 12); 
ParametersofBarrier(:,8:11) = BarrierData;
ParametersofBarrier(:,1) =  (BarrierData(:,1)+BarrierData(:,3))/2;% position of fractures
ParametersofBarrier(:,2) =  (BarrierData(:,2)+BarrierData(:,4))/2;% position of fractures
ParametersofBarrier(:,3) =  sqrt((BarrierData(:,1)-BarrierData(:,3)).^2+(BarrierData(:,2)-BarrierData(:,4)).^2);% position of fractures
ParametersofBarrier(:,4) =  atan2(BarrierData(:,4)-BarrierData(:,2),BarrierData(:,3)-BarrierData(:,1));% position of fractures


figure;surf(reshape(coordinates(:,1),Cell_N+1,Cell_M+1),reshape(coordinates(:,2),Cell_N+1,Cell_M+1),zeros(Cell_N+1,Cell_M+1));
colormap('white');axis image;axis off; ax=axis;axis(ax*1.001); view([0,90]); hold on;
plot(ParametersofFracture(:,[8,10])',ParametersofFracture(:,[9,11])',...
   'r-','LineWidth',2); 
plot(ParametersofBarrier(:,[8,10])',ParametersofBarrier(:,[9,11])',...
   'b-','LineWidth',2); hold off; 

% case (a)
%text(0-0.077,0.5,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
%text(1+0.077,0.5,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
%text(0.5,1+0.047,'p_D=4','HorizontalAlignment','center','Color','k','FontSize',16)
%text(0.5,0-0.042,'p_D=1','HorizontalAlignment','center','Color','k','FontSize',16)

% case (b)
text(0-0.077,0.5,'p_D=4','HorizontalAlignment','center','Color','k','FontSize',16)
text(1+0.077,0.5,'p_D=1','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,1+0.047,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)
text(0.5,0-0.042,'q_N=0','HorizontalAlignment','center','Color','k','FontSize',16)

