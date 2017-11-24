function aus_abgeleitet=dh_dw_PK(temp,konz)
% temp                     temperature            	K (Kelvin)
%                                                	[ ... ]
% konz                     mass fraction         	kg LiBr/kg solution 
%                                                	[0.4 ... 0.7]        
%                                                
%
% This function computes the partial derivative of the enthalpy of aqueous LiBr solution based on a funktion 
% published by Patek and Klomfar
% Input parameters are: temperature in Kelvin and mass fraction of LiBr in kg LiBr/kg solution.
% Source: J. Patek & J. Klomfar, A computationally effective formulation 
%         of the thermodynamic properties of LiBr-H2O solutions from 
%         273 to 500 K over full composition range. Published in IJR
% 
% The function returs the partial derivative of the enthalpy in kJ/kg
% The variables 'temp' and 'konz' may be vector but must be of same size

% Change log 
% 2017-11-23 	first release		M.Mittermaier



beta(1)=1/3;
beta(2)=2/3;
beta(3)=5/6;
beta(4)=21/6;
alpha(1)=-4.37196E-1;
alpha(2)=3.0344E-1;
alpha(3)=-1.29582;
alpha(4)=-1.7641E-1;
T_c=647.096;
h_c=37548.5;
m(1)=1;
m(2)=1;
m(3)=2;
m(4)=3;
m(5)=6;
m(6)=1;
m(7)=3;
m(8)=5;
m(9)=4;
m(10)=5;
m(11)=5;
m(12)=6;
m(13)=6;
m(14)=1;
m(15)=2;
m(16)=2;
m(17)=2;
m(18)=5;
m(19)=6;
m(20)=7;
m(21)=1;
m(22)=1;
m(23)=2;
m(24)=2;
m(25)=2;
m(26)=3;
m(27)=1;
m(28)=1;
m(29)=1;
m(30)=1;
n(1)=0;
n(2)=1;
n(3)=6;
n(4)=6;
n(5)=2;
n(6)=0;
n(7)=0;
n(8)=4;
n(9)=0;
n(10)=4;
n(11)=5;
n(12)=5;
n(13)=6;
n(14)=0;
n(15)=3;
n(16)=5;
n(17)=7;
n(18)=0;
n(19)=3;
n(20)=1;
n(21)=0;
n(22)=4;
n(23)=2;
n(24)=6;
n(25)=7;
n(26)=0;
n(27)=0;
n(28)=1;
n(29)=2;
n(30)=3;
t(1)=0;
t(2)=0;
t(3)=0;
t(4)=0;
t(5)=0;
t(6)=1;
t(7)=1;
t(8)=1;
t(9)=2;
t(10)=2;
t(11)=2;
t(12)=2;
t(13)=2;
t(14)=3;
t(15)=3;
t(16)=3;
t(17)=3;
t(18)=3;
t(19)=3;
t(20)=3;
t(21)=4;
t(22)=4;
t(23)=4;
t(24)=4;
t(25)=4;
t(26)=4;
t(27)=5;
t(28)=5;
t(29)=5;
t(30)=5;
a(1)=2.27431;
a(2)=-7.99511;
a(3)=3.85239E2;
a(4)=-1.6394E4;
a(5)=-4.22562E2;
a(6)=1.13314E-1;
a(7)=-8.33474;
a(8)=-1.73833E4;
a(9)=6.49763;
a(10)=3.24552E3;
a(11)=-1.34643E4;
a(12)=3.99322E4;
a(13)=-2.58877E5;
a(14)=-1.93046E-3;
a(15)=2.80616;
a(16)=-4.04479E1;
a(17)=1.45342E2;
a(18)=-2.74873;
a(19)=-4.49743E2;
a(20)=-1.21794E1;
a(21)=-5.83739E-3;
a(22)=2.3391E-1;
a(23)=3.41888E-1;
a(24)=8.85259;
a(25)=-1.78731E1;
a(26)=7.35179E-2;
a(27)=-1.7943E-4;
a(28)=1.84261E-3;
a(29)=-6.24282E-3;
a(30)=6.84765E-3;
T_0=221;
Mw=18.01528;
Ms=6.941+79.904;

if((ndims(temp)==2)&&(ndims(konz)==2)&&(ismatrix_mod(temp)==0)&&(ismatrix_mod(konz))==0)
    if (iscolumn_mod(temp))
        temp=temp';
    end
    if (iscolumn_mod(konz))
        konz=konz';
    end
    alpha=alpha';
    beta=beta';
    a=a';
    m=m';
    n=n';
    t=t';
    x = (konz/Ms)./(konz/Ms + (1-konz)/Mw);
    x_abgeleitet=((1/Ms)*(konz/Ms + (1-konz)/Mw)-(konz/Ms)*(1/Ms-1/Mw))./(konz/Ms + (1-konz)/Mw).^2;

    alpha_2=repmat(alpha,1,length(temp));
    beta_2=repmat(beta,1,length(temp));
    temp_2=repmat(temp,length(alpha),1);
    h_strich=h_c*(1+sum(alpha_2.*(1-temp_2/T_c).^beta_2));

    a_2=repmat(a,1,length(temp));
    m_2=repmat(m,1,length(temp));
    n_2=repmat(n,1,length(temp));
    t_2=repmat(t,1,length(temp));
    temp_3=repmat(temp,length(a),1);
    x_2=repmat(x,length(a),1);
    h=(1-x).*h_strich+h_c*sum(a_2.*x_2.^(m_2).*(0.4-x_2).^(n_2).*(T_c./(temp_3-T_0)).^(t_2));
    
    x_abgeleitet_2=repmat(x_abgeleitet,length(a),1);
    for k=1:30
        if(n(k)~=0)
            h_abgeleitet=-h_strich.*x_abgeleitet+h_c*sum(a_2.*(T_c./(temp_3-T_0)).^(t_2).*x_abgeleitet_2.*(m_2.*x_2.^(m_2-1).*(0.4-x_2).^(n_2)-x_2.^(m_2).*n_2.*(0.4-x_2).^(n_2-1)));
        else
            h_abgeleitet=-h_strich.*x_abgeleitet+h_c*sum(a_2.*(T_c./(temp_3-T_0)).^(t_2).*x_abgeleitet_2.*(m_2.*x_2.^(m_2-1).*(0.4-x_2).^(n_2)));
        end
    end

    M_ges=(x*Ms+(1-x)*Mw)/1000;
    M_ges_abgeleitet=(Ms*x_abgeleitet-Mw*x_abgeleitet)/1000;

    aus_abgeleitet=-(h_abgeleitet.*(1000*M_ges)-h*1000.*M_ges_abgeleitet)./(1000*M_ges).^2;
else
    warning('vectors but must be of same size!')
    aus_abgeleitet=NaN;
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
