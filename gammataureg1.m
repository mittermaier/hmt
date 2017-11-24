function gammataureg1=gammataureg1(tau, pic)
% gammataureg1(tau, pic)
%
% tau    reduced temperature     dimensionless
% pic    reduced pressure        dimensionless
%
% IAPWS water properties 1997

% Unterfunktion zu den Stoffwerten für Wasser
% Stefan Petersen 10/10/02

Stoffparameter_H2O;

tauvec(1:length(tau),1:length(nreg1))=repmat(tau-1.222,1,length(nreg1));
picvec(1:length(pic),1:length(nreg1))=repmat(7.1-pic,1,length(nreg1));

nreg1=repmat(nreg1,length(tau),1);
ireg1=repmat(ireg1,length(tau),1);
jreg1=repmat(jreg1,length(tau),1);

gammataureg1_tmp=nreg1.*picvec.^ireg1.*jreg1.*tauvec.^(jreg1-1);

gammataureg1=sum(gammataureg1_tmp')';

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
