function gamma0taureg2=gamma0taureg2(tau, pic)
% First derivative in tau of ideal-gas part of fundamental equation for region 2
% Stefan Petersen
% 10/10/02

Stoffparameter_H2O;

% tau=ones(1,9)*tau;
% gamma0taureg2=sum(n0reg2.*j0reg2.*tau.^(j0reg2-1));



tauvec(1:length(tau),1:length(n0reg2))=repmat(tau,1,length(n0reg2));

n0reg2=repmat(n0reg2,length(tau),1);
j0reg2=repmat(j0reg2,length(tau),1);


gamma0taureg2_tmp=n0reg2.*j0reg2.*tauvec.^(j0reg2-1);

gamma0taureg2=sum(gamma0taureg2_tmp')';

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
