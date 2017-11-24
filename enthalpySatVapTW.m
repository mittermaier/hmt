function  enthalpySatVapTW=enthalpySatVapTW(temp)
%  enthalpySatVapTW=enthalpySatVapTW(temp)
%
%  specific enthalpy of saturated steam as a function of temperature
%
%  enthalpySatVapTW   in kJ/kg
%  temp               in K


press = pSatW(temp);
% region 2
t2=find(temp>=273.15 & temp<=623.15);
enthalpySatVapTW(t2,1) = enthalpyreg2(temp(t2), press(t2));

% region 3
% t3=find(temp>623.15 & temp<=647.096);
%         pressure = pSatW(temperature) - 0.00001
%         density = densreg3(temperature, pressure)
%         enthalpySatVapTW = enthalpyreg3(temperature, density)

% outside range
tbad=find(temp<273.15 | temp>623.15);
enthalpySatVapTW(tbad,1) = -1;
if ~isempty(tbad)
  disp('ERROR in function enthalpySatVapTW!!! Temperature outside range. Enthalpy is set to -1');
end

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
