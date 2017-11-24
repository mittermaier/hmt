function aus=longitudgridpnt_GS_2(z_vektor,n,nu,g,u0)
% z_vektor		vector containing transversal grid points           	in m
%                                                				[ ... ]
% n			index of sought longitudinal grid point
%                                                				[ ... ]
% nu 			kinematic viskosity 					in m^2/s
%
% g 			acceleration of gravity 				in m/s^2
%										[9.81]
%u0 			initial velocity 					in m/s
%										[0.1...0.5]
%
%This function calculates longitudinale grid points. For the index n 
%the longitudinal grid point is adjusted in a way
%that the boundary layer thickness fits to the transversal grid point of the 
%given vector z_vektor at index n.
%Therefore the root of function 'Filmdicke_hilf' is sought.  

%x is the sought root. It is expected to be between 10^-19 and 1.
%10^-19 is chosen close to zero, but 0 would lead to an error in matlab. 

%The function was programmed 2012 by P.Schulze during his bachelor thesis 
%under supervision of M.Mittermaier 

% Change log 
% 2017-11-23 	first release		M.Mittermaier

aus=fzero(@(x)Filmdicke_hilf(z_vektor(1,n),nu,u0,g,x),[10^-19 1]);

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


