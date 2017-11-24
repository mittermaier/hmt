function longitudgridpnt=longitudgridpnt_GS_exp_lin(n_3,erstes_Intervall,Faktor)
% n_3                      determines which grid point is calculated            	
%                                                	
% erstes_Intervall         initial step size          	 
%                                                	
% Faktor		   base of exponential step size growth
%
% This function calculates longitudinal grid points.
% Starting from a given initial step size, the grid size grows exponentially. 
% Once the step size exceeds 100 times the initial one, all following steps stay of
% equal size 
%The function was programmed 2012 by P.Schulze during his bachelor thesis 
%under supervision of M.Mittermaier 

% Change log 
% 2017-11-23 	first release		M.Mittermaier

x = erstes_Intervall;
n=1;%counter of while loop
%in vector z the calculated coordinates are saved
z(1,1)=0;%z starts at 0
z(2,1)=x;
y=x;%y is substitute for individual values of z in subsequent loops
%Exponential part:
while (z(n+1,1)-z(n,1)<100*erstes_Intervall) %exponential growth until
%100times of the initial ostep size
    y=y+x*Faktor^n;% longitudinal step size grows 
    z(n+2,1)=y;%vector z gets value from substitute y
    n=n+1;%counter
end
%equidistant part:
for m=n+1:n_3-1
    y=y+z(n+1,1)-z(n,1);%final exponentially generated step size
    z(m+1,1)=y;%vector z gets value from substitute y
end
%The output is the value at step n_3 (n_3 is transferred as input)
longitudgridpnt = z(n_3,1);

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
