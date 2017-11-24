function enthalpyreg1=enthalpyreg1(temp, press)
% specific enthalpy in region 1
% Spezifische Enthalpie für Wasser im flüssigen Zustand
% enthalpyreg1  in kJ/kg
% temp          in K
% press         in bar
%
% erstellt von Stefan Petersn
% 10/10/02
%

Stoffparameter_H2O;
tau=1386 ./ temp;
pic=0.1 .* press ./ 16.53;
enthalpyreg1 = 0.001 .* rgas_water .* temp .* tau .* gammataureg1(tau, pic);
