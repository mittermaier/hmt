function enthalpySatLiqTW=enthalpySatLiqTW(temp)
% enthalpySatLiqTW=enthalpySatLiqTW(temp)
% specific enthalpy of saturated liquid water as a function of temperature
%
% enthalpySatLiqTW    in kJ/kg
% temp                in K
%

press = pSatW(temp);
% region 1
t1=find(temp>=273.15 & temp<=623.15);
enthalpySatLiqTW(t1,1) = enthalpyreg1(temp(t1), press(t1));

% region X

% outside range
tbad=find(temp<273.15 | temp>623.15);
enthalpySatLiqTW(tbad,1) = -1;
if ~isempty(tbad)
  disp('ERROR in function enthalpySatLiqTW!!! Temperature outside range. Enthalpy is set to -1');
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
