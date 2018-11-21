function [G,V]=gradvecglmalphaup_noZonal(TH,L,rnew,rold)

  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);
  
  % The spherical harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;

  % Find row indices into G belonging to the orders
  [EM,EL,mz,blkm]=addmout(maxL);
  
  % Find increasing column index; that's how many belong to this order
  % alpha=cumsum([1 L+1 gamini(L:-1:1,2)]);
  % The middle bit is twice for every nonzero order missing
  % alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  %   		gamini(L(2-lp)-bp*(L(1)-1),bp*2*(L(1)-1)) ...
  %   		gamini(L(2-lp)-bp*(L(1)-1):-1:1,2)]);
  % This should be the same for L and [0 L]
  alpha=cumsum([1 L(2-lp)-bp*L(1)+1 ...
  		gamini(L(2-lp)-bp*(L(1)-1),bp*2*L(1)) ...
  		gamini(L(2-lp)-bp*L(1):-1:1,2)]);
  
  % For AXISYMMETRIC REGIONS
  % Initialize matrices
  G=zeros((maxL+1)^2,ldim);%repmat(0,(maxL+1)^2,ldim);
  V=zeros(1,ldim);
  disp('Calculating in parallel mode')
  try
    parpool
  end
  
  parfor mm=1:maxL+1       
    m=mm-1;
    %[E,Vpp,Np,th,Cp]=gradvecsdwcapup(TH,L,m,0,-1,[],[],rnew,rold);
    [~,Vpp,~,~,Cp]=gradvecsdwcapup(TH,L,m,0,-1,[],[],rnew,rold);
    Vp{mm}=Vpp;
    C{mm}=Cp;
  end 
  % Distribute this at the right point in the huge matrix
  % Put the m=0 as nan
  m=0;
  G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=nan(size(C{m+1}));
  V(alpha(2*m+1):alpha(2*m+2)-1)=nan(size(Vp{m+1}));
  
  for m=1:maxL%m=0:maxL
    if m>0
      % Here you supply the negative orders
      G(EM==-m,alpha(2*m):alpha(2*m+1)-1)=C{m+1};
      V(alpha(2*m):alpha(2*m+1)-1)=Vp{m+1};
      %MTAP(alpha(2*m):alpha(2*m+1)-1)=-m;
      % It's all neatly ordered here, downgoing within every order
      %IMTAP(alpha(2*m):alpha(2*m+1)-1)=1:length(Vp{m+1});
    end
    % Duplicate for the positive order in case the region is axisymmetric
    G(EM==m,alpha(2*m+1):alpha(2*m+2)-1)=C{m+1};
    V(alpha(2*m+1):alpha(2*m+2)-1)=Vp{m+1};
    %MTAP(alpha(2*m+1):alpha(2*m+2)-1)=m;
    % It's all neatly ordered here, downgoing within every order
    %IMTAP(alpha(2*m+1):alpha(2*m+2)-1)=1:length(Vp{m+1});
  end

  % Remove those rows with m=0 and the first maxL columns
  G(mz,:)=[];
  G(:,1:maxL+1)=[];
  V(1:maxL+1)=[];
  
  [V,isrt]=sort(V,'descend');
  G=G(:,isrt);
  

