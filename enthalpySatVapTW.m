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
