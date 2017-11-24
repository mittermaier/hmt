%  This M-file contains the coefficients for calculating the thermodynamical properties 
%  of LiBr-H2O solutions from 273 to 500 K over full composition range
%
%  Source: J. Patek & J. Klomfar, A computationally effective formulation 
%          of the thermodynamic properties of LiBr-H2O solutions from 
%          273 to 500 K over full composition range. To be published in IJR


% Critical point of pure water
TCritW          = 647.096;          % K
pCritW          = 22.064*10^6;      % Pa
densCritW       = 322 ;             % kg/m³
densCritWmol    = 17873.727 ;       % mol*m³
cpCritWmol      = 76.0226 ;         % J/(mol-K)
enthalpyCritWmol= 37548.5;          % J/mol
entropyCritWmol = 79.3933;          % J/(mol-K)

% Triple point of pure water
TTripW = 273.16;        % K
pTripW = 611.657;       % Pa
densTripW = 999.789 ;   % kg/m³

% Molar mass
M_W    = 0.018015268; % kg/mol
M_LiBr = 0.08685;     % kg/mol

% Regression coefficients

% Table 4: Pressure
mTab4=[3 4 4 8 1 1 4 6];
nTab4=[0 5 6 3 0 2 6 0];
tTab4=[0 0 0 0 1 1 1 1];
aTab4=[-241.303 19175000 -175521000 32543000 392.571 -2126.26 185127000 1912.16];

% Table 5: Density
mTab5=[1 1];
tTab5=[0 6];
aTab5=[1.746 4.709];

% Table 6: cp
mTab6=[2 3 3 3 3 2 1 1];
nTab6=[0 0 1 2 3 0 3 2];
tTab6=[0 0 0 0 0 2 3 4];
aTab6=[-14.2094 40.4943 111.135 229.98 1345.26 -0.014101 0.0124977 -0.000683209];

% Table 6: Enthalpy
mTab7=[1 1 2 3 6 1 3 5 4 5 5 6 6 1 2 2 2 5 6 7 1 1 2 2 2 3 1 1 1 1];
nTab7=[0 1 6 6 2 0 0 4 0 4 5 5 6 0 3 5 7 0 3 1 0 4 2 6 7 0 0 1 2 3];
tTab7=[0 0 0 0 0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5];
aTab7=[2.27431 -7.99511 385.239 -16394 -422.562 0.113314 -8.33474 -17383.3 ...
       6.49763 3245.52 -13464.3 39932.2 -258877 -0.00193046 2.80616 -40.4479 ...
       145.342 -2.74873 -449.743 -12.1794 -0.00583739 0.23391 0.341888 8.85259 ...
       -17.8731 0.0735179 -0.00017943 0.00184261 -0.00624282 0.00684765];

% Table 8: Entropy
mTab8=[1 1 2 3 6 1 3 5 1 2 2 4 5 5 6 6 1 3 5 7 1 1 1 2 3 1 1 1 1];
nTab8=[0 1 6 6 2 0 0 4 0 0 4 0 4 5 2 5 0 4 0 1 0 2 4 7 1 0 1 2 3];
tTab8=[0 0 0 0 0 1 1 1 2 2 2 2 2 2 2 2 3 3 3 3 4 4 4 4 4 5 5 5 5];
aTab8=[1.53091 -4.52564 698.302 -21666.4 -1475.33 0.0847012 -6.59523 ...
        -29533.1 0.00956314 -0.188679 9.31752 5.78104 13893.1 -17176.2 ...
        415.108 -55564.7 -0.00423409 30.5242 -1.6762 14.8283 0.00303055 ...
        -0.040181 0.149252 2.5924 -0.177421 -0.000069965 0.000605007 ...
        -0.00165228 0.00122966];


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
