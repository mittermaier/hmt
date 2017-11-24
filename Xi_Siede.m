function aus=Xi_Siede(T,p)
% T			temperature         		C (celcius)
%                                                	
% p                    	pressure	         	Pa 
%                                                	
% Xi_Siede		saturation mass fraction	kg LiBr/kg solution 
%
% This function calculates the saturation mass fraction (mass fraction at thermodynamic equilibrium) of aqueous LiBr
% solution. Input parameters are: mass fraction of LiBr and the pressure 
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

epsilon=1e-2;
Startwert=0.01; %initial guess value; 0.0 not defined!
x0=Startwert;
x1=Startwert;

if T_Siede(Startwert,p)-T > 0
   while x0==Startwert
       if T_Siede(x1+epsilon,p)-T > 0
           x1=x1+epsilon;
       end
       if T_Siede(x1+epsilon,p)-T < 0
           x0=x1+epsilon;
       end
   end
end

if T_Siede(Startwert,p)-T < 0
    while x1==Startwert
        if T_Siede(x0+epsilon,p)-T < 0
            x0=x0+epsilon;
        end
        if T_Siede(x0+epsilon,p)-T > 0
            x1=x0+epsilon;
        end
    end
end

xalt=x0;
xaltalt=x1;
xneu=x1; 
durchlauf=0;
while abs(xneu-xalt) > 1e-7
   durchlauf=durchlauf+1;
   if durchlauf > 1
       xaltalt=xalt;
       xalt=xneu;
   end
   if abs((T_Siede(xalt,p)-T)-(T_Siede(xaltalt,p)-T))>0
       xneu=xalt-(xalt-xaltalt)/((T_Siede(xalt,p)-T)-(T_Siede(xaltalt,p)-T))*(T_Siede(xalt,p)-T);
   else
       break;
   end
end
aus=xneu;


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
