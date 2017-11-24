function aus=transversalgridpnt(y_Anzahl,delta)

% y_Anzahl			number of nodes or grid points representing the film in 
%                               transversal direction.
%
% delta				film thickness 			in meter         
%                                                
%
% This function generates transversal grid coordinates and returns them as a vector. 
% Input parameters are: number of gridpoints 'y_Anzahl' and the film thickness 'delta' 
% they should represent.
%
% The grid points (nodes) are NOT equally spaced. The grid spaces are small towards the wall and
% towards the interface. The maximum space between the nodes is in the middle of the film thickness. 
% The distribution follows a bell curve (normal distribution) that is scaled to the film thickness. 

% Change log 
% 2017-11-23 	first release		M.Mittermaier


sig=delta/4; 	% is an arbitratry value should be between 4 and 8 
		% with regular film thicknesses of 0.00025m the order of magnitude of the spaced is  10^-10 to 10^-7
              	% the smaller variable 'sig', the larger the difference in grid space 
               

y_koord=0:delta/y_Anzahl:delta; %vector with coordinates

funk_wert= 1/(sig*sqrt(2*pi)) .* exp(-0.5*((y_koord - delta/2)./sig).^2);

funk_wert = (funk_wert(1:y_Anzahl)+funk_wert(2:y_Anzahl+1))*0.5;

skal=delta/sum(funk_wert,2); %scaling the bell curve

aus=skal.*funk_wert;
