function aus=d_vec(T,Xi)
% T                        temperature            	C (Celcius)
%                                                	[ ... ]
% Xi                       mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]
% D			   mass diffusivity     	in m^2/s
%
% This function calculates the mass diffusivity of aqueous LiBr
% solution. Input parameters are: temperature and mass fraction of LiBr
%
% The calculation is based on the PHD thesis p.189: 

% Kim, Kwang J
% Heat and mass transfer enhancement in absorption cooling 
% Arizona State University
% 1992

% and personal communication between Kwang J. Kim and Felix Ziegler


% Change log 
% 2017-11-23 	Made suitable for vectors, first release		M.Mittermaier

m=1000*Xi/86.845./(1-Xi);
D25=(1.3528+0.19881*m-0.036382*m.^2+0.0020299*m.^3-0.000039375*m.^4)*1e-9;
D=D25.*(T+273.15)./298.15.*eta_LiBr_vec(25,Xi)./eta_LiBr_vec(T,Xi);
aus=D;


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
