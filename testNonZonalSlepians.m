function testNonZonalSlepians()

  index=10;
  cmp=3;
  Lmax=20;
  dom=10;

  rplanet=2440;
  rsat=2540;
  
  % Calculate Slepian coefficients
  [H,V]=gradvecglmalphaup_noZonal(dom,Lmax,rsat,rplanet);


  

%%% Evaluate at altitude %%%
  % Insert the zonal coefficients as zeros to make evaluation easier
  Hz=insertZones(H,Lmax);
  % Upward continuation
  Hz=vecupderivative(Hz,rsat,rplanet,Lmax);
  % Transform the index-th best
  elmcosi=coef2lmcosi(Hz(:,index),1);
  % Evaluate the Slepian functions
  [r,lon,lat]=elm2xyz(elmcosi,1,[0 89 360 -89]);

  
%%% Now evaluate with rGvec %%%
  [LON,LAT]=meshgrid(lon,lat);
  rG=rGvec_noZonal(H(:,index),(90-LAT(:))*pi/180,LON(:)*pi/180,...
                   rsat*ones(size(LON(:))),rplanet,Lmax,1);
  

  subplot(1,3,1)
  plotplm(r{cmp}/sqrt(4*pi),lon*pi/180,lat*pi/180,2,1)
  view(90,90)
  kelicol(1)
  cval=max(abs(caxis));
  caxis([-1,1]*cval)

  subplot(1,3,2)
  %keyboard
  rGcomp=reshape(rG((cmp-1)*length(LON(:))+1:cmp*length(LON(:))),length(lat),length(lon));
  plotplm(rGcomp,lon*pi/180,lat*pi/180,2,1)
  kelicol(1)
  caxis([-1,1]*cval)
  view(90,90)

  subplot(1,3,3)
  diff=rGcomp-r{cmp}/sqrt(4*pi);
  plotplm(diff,lon*pi/180,lat*pi/180,2,1)
  kelicol(1)
  caxis([-1,1]*cval)
  view(90,90)
  max(abs(diff(:)))
