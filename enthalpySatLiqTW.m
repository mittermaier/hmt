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



