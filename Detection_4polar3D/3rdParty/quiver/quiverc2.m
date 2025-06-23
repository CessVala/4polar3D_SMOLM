function quiverc2(x,y,u,v,anyso,imagegray,varargin)
global zoomxy
% This function is inspired from vfield; I tried to overcome some problems
% that I found using the other one and to improve the visulaization of the
% vectors of the field. Futhermore the visualization is CFD style.

% The visualization is very quick, also with large vector fields.

% The function reads matrix x,y positions and u,v vectors components. All
% input must have same dimensions

% The dimension of plotted vector is proportional to vector length, and the
% color is same.

% You can specify two options
% 'scale' and the value of che scale factor enlarging vector visualization
% 'lungh' and you have all vector of normalized length 


lungh=0.1*zoomxy; %normal 0.1 ... 0.5 for single one dna ...0.2 for THT and amyloid, 
scale=1;

ampl=0.1;
u=ampl*u;v=ampl*v;

for i=1:length(varargin)
  vin=varargin{i};
  if isequal(vin,'lungh')
    lungh=varargin{i+1};
  end
  if isequal(vin,'scale')
    scale=varargin{i+1};
  end
end

 speed=(u.^2+v.^2).^0.5;

[theta,r]=cart2pol(u,v);

% pause
apertura_freccia=30*pi/180;

lunghezza_freccia=8/10;

largh_freccia=0.05;

delta_primo=atan((1-lunghezza_freccia)/lunghezza_freccia*tan(apertura_freccia));
delta_freccia=atan(largh_freccia/lunghezza_freccia);

if lungh==0
    r=r*scale;
    r_primo=(lunghezza_freccia*speed/cos(delta_primo))*scale;
    r_primo_freccia=lunghezza_freccia*speed/cos(delta_freccia)*scale;
    r_largh_freccia=largh_freccia*speed*scale;
else
    r=lungh*scale*ones(size(speed));
    r_primo=(lunghezza_freccia/cos(delta_primo))*lungh*scale*ones(size(speed));
    r_primo_freccia=lunghezza_freccia/cos(delta_freccia)*lungh*scale*ones(size(speed));
    r_largh_freccia=largh_freccia*lungh*scale*ones(size(speed));
end

% [X1,Y1]=pol2cart((theta+pi/2),r_largh_freccia);
% [X2,Y2]=pol2cart((theta+delta_freccia),r_primo_freccia);
% [X3,Y3]=pol2cart((theta+delta_primo),r_primo);
% [X4,Y4]=pol2cart(theta,r);
% [X5,Y5]=pol2cart((theta-delta_primo),r_primo);
% [X6,Y6]=pol2cart((theta-delta_freccia),r_primo_freccia);
% [X7,Y7]=pol2cart((theta-pi/2),r_largh_freccia);

newr=atan(r_largh_freccia./r);
newl=sqrt(r_largh_freccia.^2+r.^2);

newr2=atan(r_largh_freccia./(r./2));
newl2=sqrt(r_largh_freccia.^2+(r./2).^2);

[X1,Y1]=pol2cart((theta+pi/2),r_largh_freccia);
[X2,Y2]=pol2cart((theta+newr2),newl2);

[X3,Y3]=pol2cart((theta+newr),newl);

[X4,Y4]=pol2cart(theta,r);

[X5,Y5]=pol2cart((theta-newr),newl);

[X6,Y6]=pol2cart((theta-newr2),newl2);
[X7,Y7]=pol2cart((theta-pi/2),r_largh_freccia);


Vert=[[x y]+[X1 Y1];
    [x y]+[X2 Y2];
    [x y]+[X3 Y3];
    [x y]+[X4 Y4];
    [x y]+[X5 Y5];
    [x y]+[X6 Y6];
    [x y]+[X7 Y7];
];

%Move to center

dx=-(r./2).*cos(theta);
dy=-(r./2).*sin(theta);
Vert=[[x+dx y+dy]+[X1 Y1];
    [x+dx y+dy]+[X2 Y2];
    [x+dx y+dy]+[X3 Y3];
    [x+dx y+dy]+[X4 Y4];
    [x+dx y+dy]+[X5 Y5];
    [x+dx y+dy]+[X6 Y6];
    [x+dx y+dy]+[X7 Y7];
];



VertX=Vert(:,1:size(Vert,2)/2);
VertY=Vert(:,size(Vert,2)/2+1:end);

VertX=reshape(VertX',[],7);
VertY=reshape(VertY',[],7);

Color=reshape((anyso)',[],1);

% figure(1);
hold on
 set(gcf,'color','w');
 set(gca,'color','k');
 imagesc(imagegray)
%  alpha(0.9)
 hold on
axis image
% axis equal
patch(VertX',VertY',Color','edgecolor','none');
% return
set(gca,'box','on')
colorbar     
set(colorbar,'Xcolor','k','Ycolor','k','Zcolor','k')
caxis([-0.6 0.6])





