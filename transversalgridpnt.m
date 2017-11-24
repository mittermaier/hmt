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

%Copyright 2017 Martin Mittermaier
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
%    der GNU General Public License, wie von der Free Software Foundation,
%    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
%    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
%
%    Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
%    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
%    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
%    Siehe die GNU General Public License für weitere Details.
%
%    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
%    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
