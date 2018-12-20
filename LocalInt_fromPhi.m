function coef=LocalInt_fromPhi(data,rad,cola,lon,dom,Lmax,J,rplanet)
% data      vectorial data:
%              given as data{1}=rad component, data{2}=colat component,
%              data{3}=lon component, 
%           Components can be ommitted: data{1}=[] for example.
% rad       radial position of satellite (planet radius + altitude)
% cola,lon  colatitude/longitude positions of the data values (both as
%           column values), 0<=cola<=pi; 0<=lon<=2pi
% dom       integer (polar cap opening angle [degrees])
% Lmax      Maximum spherical harmonic degree for the internal source part
% J         How many Slepian functions should be used to calculate the
%           solution? More means more sensitive to noise but higher spatial
%           resolution
% rplanet   Planet radius to which the inner source solution should be 
%           calculated
%
% Last modified, 12/19/2018 plattner-at-alumni.ethz.ch
  
  useit(1)=~isempty(data{1});
  useit(2)=~isempty(data{2});
  useit(3)=~isempty(data{3});
  data=[data{1};data{2};data{3}];
  
  avgsat=mean(rad);
  [H,V]=gradvecglmalphaup_noZonal(dom,Lmax,avgsat,rplanet);

  % mz is where the m=0 sit.
  [~,~,mz,~]=addmout(max(Lmax));

  % Now evaluate
  MloadJ=rGvec_noZonal(H(:,1:J),cola,lon,rad,rplanet,max(Lmax),1);
  MloadJ=MloadJ(:,repelem(useit,length(rad)));

  M=MloadJ(1:J,:);

  MM=M*M';
  Md=M*data;
  slepcoef=MM\Md;

  % And now iteratively reweighted residual calculation
  [slepcoef,dataweights]=itweighres(M,data,slepcoef,10);
  
  coef=H(:,1:J)*slepcoef;  
  
  % Fill in zeros for zones
  coef=insertZones(coef,max(Lmax));

  % Coefs in addmout. Transform to addmon
  coef=out2on(coef,max(Lmax));
