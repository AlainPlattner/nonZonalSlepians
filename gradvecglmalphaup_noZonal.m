function [G,V]=gradvecglmalphaup_noZonal(TH,L,rnew,rold)
  % [G,V]=gradvecglmalphaup(TH,L,rnew,rold,srt,anti)
  % 
  % Construction of gradient vector Slepian functions
  %
  % INPUT:
  % 
  % TH    Region name (e.g. 'africa' etc),
  %       OR opening angle of spherical cap,
  %       OR Angles of two spherical caps and we want the ring between 
  %       them [TH1 TH2]
  %       OR [lon lat] an ordered list defining a closed curve [degrees]
  %       OR several regions to add up/subtract: 
  %          struct with
  %          TH.name    for name of the combined region
  %          TH.parts   for the cell array of names of the parts, 
  %                     or cap opening angles
  %                     or [cap,lon,colat] for rotated caps
  %          TH.sign    for adding or subtracting 
  %          Example: TH.parts{1}='namerica'; TH.parts{2}='samerica';
  %                   TH.sign=[1,1]; TH.name='americas';
  %                   TH.name='weirdRing'
  %                   TH.parts{1}=30; TH.parts{2}=[5,5,10]; TH.sign=[1,-1]
  %                   subtracts the ring of cTH=5, clon=5, ccola=10 from the
  %                   larger polar cap
  % L     Bandwidth (maximum angular degree), or passband (two degrees)
  % rnew  Satellite altitude
  % rold  planet radius
  %
  % OUTPUT:
  %
  % G     Matrix containing vector Slepian function coefficients for the Elm
  %       vector spherical harmonics IN ADDMOUT FORMAT
  % V     Eigenvalues (conditioning values) UNSORTED FOR REGULAR REGIONS
  %
  % Spherical cap part: plattner-at-alumni.ethz.ch, 5/23/2019
  % Last modified: plattner-at-alumni.ethz.ch, 1/22/2024
  
  % Figure out if it's lowpass or bandpass
  lp=length(L)==1;
  bp=length(L)==2;
  maxL=max(L);

  defval('anti',false);
  
  % The spherical harmonic dimension
  ldim=(L(2-lp)+1)^2-bp*L(1)^2;


% First check if already calculated
  if ~isstr(TH) && ~isstruct(TH) && length(TH)==1 % POLAR CAPS
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphaup-nz-%g-%i-%g-%g.mat',TH,L,rnew,rold));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphablup-nz-%g-%i-%i-%g-%g.mat',...
             TH,L(1),L(2),rnew,rold));
    else
      error('The degree range is either one or two numbers')       
    end
    
  elseif ~isstr(TH) && ~isstruct(TH) && length(TH)==2 % Ring between polar cap angles
    defval('sord',1) % SINGLE OR DOUBLE CAP
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphaup-nz-%g_%g-%i-%g-%g.mat',max(TH),...
             min(TH),L,rnew,rold));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphablup-nz-%g_%g-%i-%i-%g-%g.mat',...
			     max(TH),min(TH),L(1),L(2),rnew,rold));
    else
      error('The degree range is either one or two numbers')       
    end

    % Initialize ordering matrices
    %MTAP=repmat(0,1,ldim);
    %IMTAP=repmat(0,1,ldim);        
  else % GEOGRAPHICAL REGIONS and XY REGIONS
    defval('sord',10) % SPLINING SMOOTHNESS
    % We'll put in a Shannon number based on the area only, not based on
    % an actual sum of the eigenvalues
    defval('J',ldim)
    % Note the next line, though we can change our minds
    %defval('J',ldim*spharea(TH))
    if isstr(TH) % Geographic (keep the string)
      h=TH;
    elseif isstruct(TH)
      h=TH.name;
    else % Coordinates (make a hash)
      try
	h=hash(TH,'sha1');
      catch
  h=builtin('hash','sha1',TH);
      end
    end
    if lp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphaup-nz-%s-%i-%i-%g-%g-%i.mat',...
             h,L,J,rnew,rold,anti));
    elseif bp
      fname=fullfile(getenv('IFILES'),'GRADVECGLMALPHAUP-NZ',...
		     sprintf('gradvecglmalphablup-nz-%s-%i-%i-%i-%g-%g-%i.mat',...
             h,L(1),L(2),J,rnew,rold,anti));
    else
     error('The degree range is either one or two numbers')       
    end
    defval('GM2AL',NaN) % If not, calculate order per taper
    defval('MTAP',NaN) % If not, calculate order per taper
    defval('IMTAP',NaN) % And rank ordering within that taper
    defval('xver',0) % For excessive verification of the geographical case
  end
  
if exist(fname,'file')==2
  load(fname)
  disp(sprintf('Loading %s',fname))
else
  
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


  % For GEOGRAPHICAL REGIONS or XY REGIONS
  if isstr(TH) || isstruct(TH) || length(TH)>2
  % Initialize matrices
  %G=zeros((maxL+1)^2,ldim);
  %V=zeros(1,ldim);  
    if bp
      error('Bandpass geographical tapers are not ready yet')
    end

        if isstruct(TH)
      % Several named regions. We will add them up.
      Klmlmp=zeros((L+1)^2,(L+1)^2);
      for reg=1:length(TH.parts)
          % If the subregion is a named region
          if ischar(TH.parts{reg})
             Kreg=kernelepup(L,TH.parts{reg},rnew,rold);
          else
             % If the subregion is a polar cap
             if length(TH.parts{reg})==1
                 % North-polar cap
                 Kreg=kernelepupcap(L,TH.parts{reg},rnew,rold);
             else
                 % Cap that needs to be rotated
                 cTH=TH.parts{reg}(1);
                 rotlon=TH.parts{reg}(2);
                 rotcola=TH.parts{reg}(3);
                 Kreg=kernelepupcap(L,cTH,rnew,rold,[rotlon,rotcola]);
             end
          end
        Klmlmp=Klmlmp + TH.sign(reg)*Kreg;
        
      end
      
    else 
        % If it's not a struct, it's either a string or a list of
        % coordinates. In both cases, kernelepup takes care of it.
        Klmlmp=kernelepup(L,TH,rnew,rold,[],[],[],anti);
    end


%%%%% Need to remove the zonal entries from K
%%% Klmlmp is in addmon.
        [~,~,~,~,~,~,bigm]=addmon(L);
        zind = bigm==0;
        Klmlmp(zind,:)=[];
        Klmlmp(:,zind)=[];

        % Calculates the eigenfunctions/values for this localization problem
        [Gnz,V]=eig(Klmlmp);
        [V,isrt]=sort(sum(real(V),1));
        V=fliplr(V);
        Gnz=Gnz(:,fliplr(isrt));

        %%% For switch to addmout: put in zero terms, switch, remove them...
%%%%%% Put zeros in for the zonal terms
        mnzind = bigm~=0;       
        %%% To do that: Make full matrix of zeros, then put the non-zonal
        %%% terms back in
        G = zeros(length(mnzind),sum(mnzind));
        G(mnzind,:)=Gnz;
        %%% This is taken care of by the insertZones.m function
    
    
    [a,b,c,d,e,f,ems,els,R1,R2]=addmon(L);
    % This indexes the orders of G back as 0 -101 -2-1012 etc
    G=G(R1,:);

%%% Remove nonzonal zeros to complete switch to addmout
    %% But need the correct location for addmout zeros
    zindout = zind(R1);
    G(zindout,:)=[];
    % Check indexing
%     difer(els(R1)-EL,[],[],mesg)
%     difer(ems(R1)-EM,[],[],mesg)
    
    % Calculate Shannon number and compare this with the theory
    %N=sum(V);
    %G=G(:,1:J);
    %V=V(1:J);


    try
      % If you are running Matlab
    	save(fname,'G','V','-v7.3')
    catch
      % If you are running octave
    	save(fname,'G','V') 
    end 



  else
  
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
  

  try
    % If you are running Matlab
    save(fname,'G','V','-v7.3')
  catch
    save(fname,'G','V')  
  end 
        
  %save(fname,'G','V','EL','EM')
  end
end
