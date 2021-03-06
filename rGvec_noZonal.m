function rG=rGvec_noZonal(coefs,theta,phi,rad,rplanet,Lmax,onorout)
% rG=rGvec_noZonal(coefs,theta,phi,rad,rplanet,Lmax,onorout)
%
% Evaluates any gradient functions (e.g. gradient vector Slepian functions)
% given by the provided coefficients of Elm at the provided points. 
% This function is quite efficient as it does not assemble the entire Elm
% matrix but iterates through the degrees to save time and memory  
% 
% INPUT:
%
% coefs     (L+1)^2 -(L+1) x J matrix whose columns are
%           the (Slepian) coefficients 
% theta     colatitude of the locations in radians [0<=theta<=pi]
% phi       longitude of the locations in radians [0<=phi<=2*pi]
% rad       radial location of the data points 
% rplanet   radius for which the coefficients are defined
% Lmax      maximum spherical-harmonic degree 
% onorout   are the columns of "coefs" in ADDMON (0) or ADDMOUT (1) format?
%           default: 0
%
% OUTPUT:
%
% r         Matrix of evaluated Slepian functions, size J x 3*npoints
%           Order: radial, colatitudinal, longitudinal
%
% See also elm, rGscal
%
% Last modified by plattner-at-alumni.ethz.ch, 12/19/2018

defval('onorout',0)

% These are the m positions for the non-zonal vector
coefpos=cumsum((0:Lmax)*2);


% We need the coefficients in addmout for this to work efficiently. If they
% are in addmon, transform them to addmout
if ~onorout
   % [~,~,~,~,~,~,~,~,rinm]=addmon(Lmax);
  % coefs=coefs(rinm,:);
  error('Need coefficients in ADDMON and without Zonal parts')
end   

rG=zeros(size(coefs,2),3*length(theta));

% Make sure phi and rad are a row vectors
phi=phi(:)';
rad=rad(:)';

% Also, include the phase shift
phi=phi+pi;

divsinvals=1./sin(theta(:)');

% % Start with L=0 because it is special. for L=0 Elm is equal to Ylm:
% L=0;
% m=0;
% % Calculating the Xlm
% X=xlm(L,abs(m),theta);
% % HERE IS THE UPWARD CONTINUATION PART:
% % Multiply each Xlm with the corresponding r factor:
% X=X.*( (-L-1)/rplanet*(rad/rplanet).^(-L-2) );        
% % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
% % The negative m is the cosine
% P=diag(sqrt(2-(m(:)==0)))*...
%      cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))); 
% rG(:,1:length(theta))=nan(size(coefs(L^2+1:(L+1)^2,:)'*(X.*P))); 
 

% Iterate over the other Ls
for L=1:Lmax
    % Setting the ms
    m=-L:L;    
    % Calculating the Xlm and the dXlm
    [X,dX]=xdxlm(L,abs(m),theta);
    % Make the longitudinal phase: ones, sines or cosines, sqrt(2) or not 
    % The negative m is the cosine
    P=diag(sqrt(2-(m(:)==0)))*...
         cos(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi)));
    dP=diag(sqrt(2-(m(:)==0)))*...
	   diag(m)*(-sin(m(:)*phi-pi/2*(m(:)>0)*ones(size(phi))));
    % radial factor:
    Rfac=repmat( 1/rplanet*(rad/rplanet).^(-L-2) ,length(m),1);
    % Remove the Zonal parts
    Rfac(L+1,:)=[];
    X(L+1,:)=[];
    dX(L+1,:)=[];
    P(L+1,:)=[];
    dP(L+1,:)=[];

    % Which coefficients
    cofs = coefs(coefpos(L)+1:coefpos(L+1),:);
    
    % Formula for renormalized E 
    % (see Plattner and Simons (2015) Handbook eq. 22):
    % Erad   = (-l-1)*X*P
    % Etheta = dX*P
    % Ephi   = 1/sin(theta)*X*dP
    % And of course don't forget the radial factors.
    % Do the same thing as fin rGscal: sum up in each degree:
    % Radial component

    rG(:,1:length(theta)) =  rG(:,1:length(theta)) + ...
        cofs'*( (-L-1)*Rfac.*X.*P );
    % Colatitudinal component
    rG(:,length(theta)+1:2*length(theta)) = ...
        rG(:,length(theta)+1:2*length(theta)) + ...
        cofs'*( Rfac.*dX.*P );
    % Longitudinal component
    rG(:,2*length(theta)+1:end) = ...
        rG(:,2*length(theta)+1:end) + ... 
         cofs'* ...
         ( Rfac.*repmat(divsinvals,length(m)-1,1).*X.*dP );


end

    

