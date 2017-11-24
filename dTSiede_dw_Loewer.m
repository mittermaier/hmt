function aus=dTSiede_dw_Loewer(var1xi, var2p)
%
% var1xi                   mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% var2p 		   pressure 			Pa
%							[ ... ]
% dTSiede_dw		   derivative of saturation 
%			   temperature with respect to mass fraction
%
% This function calculates the derivative of saturation temperature with respect to mass fraction of aqueous LiBr
% solution. Input parameters are: mass fraction of LiBr and pressure 
%
% The calculation is based on data from the PHD thesis: 
% Löwer, Harald 
% Thermodynamische und physikalische Eigenschaften der wässrigen Lithiumbromid-Lösung Karlsruhe
% 1960

% The regression alaysis was conducted in 2009 during a PHD study by:
% Wohlfeil, Arnold
% Wärme- und Stoffübertragung bei der Absorption an Rieselfilmen in Absorptionskälteanlagen
% Technische Universität Berlin
% 2009

% Change log 
% 2017-11-23 	first release		M.Mittermaier

a0=-0.00470858;
a1=-0.001276757;
a2=0.000145597;
a3=0.000428261;
a4=0.000948526;
a5=3.47501E-06;
a6=-0.000495401;
a7=-5.44472E-05;
a8=0.000110477;
a9=0.004915398;
a10=-7.21234E-08;
a11=-0.00058121;
a12=-2.23738E-05;
a13=2.39788E-06;
a14=-6.64049E-06;
a15=4.26683E-06;

lnp=ln(var2p);
hilf=a0+a1*var1xi+a2*lnp+a3*var1xi.*lnp+a4*var1xi.^2+a5*lnp.^2+a6*var1xi.^2.*lnp+a7*var1xi.*lnp.^2+a8*var1xi.^2.*lnp.^2+a9*var1xi.^3+a10*lnp.^3+a11*var1xi.^3.*lnp+a12*var1xi.^3.*lnp.^2+a13*var1xi.*lnp.^3+a14*var1xi.^2.*lnp.^3+a15*var1xi.^3.*lnp.^3;
hilf_strich=a1+a3*lnp+2*a4*var1xi+2*a6*var1xi.*lnp+a7*lnp.^2+2*a8*var1xi.*lnp.^2+3*a9*var1xi.^2+3*a11*var1xi.^2.*lnp+3*a12*var1xi.^2.*lnp.^2+a13*lnp.^3+2*a14*var1xi.*lnp.^3+3*a15*var1xi.^2.*lnp.^3;
aus=-hilf_strich/(hilf^2);

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
