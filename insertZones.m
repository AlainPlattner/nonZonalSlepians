function Gfull=insertZones(G,Lmax)

  warning('G must be in ADDMOUT')
  
  Gfull=zeros((Lmax+1)^2,size(G,2));

  % mz is where the m=0 sit
  [~,~,mz]=addmout(Lmax);
  % We won't do L=0
  mz=mz(2:end);

  % coefpos are the positions for the 
  coefpos=cumsum((0:Lmax)*2);

  for L=1:Lmax
    % The negative ones
    Gfull(mz(L)-L:mz(L)-1,:)=G(coefpos(L)+1:coefpos(L)+L,:);
    % The positive ones
    Gfull(mz(L)+1:mz(L)+L,:)=G(coefpos(L+1)-L+1:coefpos(L+1),:);
  end
