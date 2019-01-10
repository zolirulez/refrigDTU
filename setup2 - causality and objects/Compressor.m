classdef Compressor   < matlab.mixin.Copyable
   properties
      Volume
      IsentropicEfficiency
   end
   methods
      function Dm = flow(comp,frequency,dInlet)
         Dm = 2*pi*frequency*dInlet*comp.Volume;
      end
      function hOutlet = enthalpy(comp,pInlet,pOutlet,hInlet)
         s = CoolProp.PropsSI('S','H',hInlet,'P',pInlet,'CO2');
         hOutletIdeal = CoolProp.PropsSI('H','S',s,'P',pOutlet,'CO2');
         hOutlet = (hOutletIdeal - hInlet)/comp.IsentropicEfficiency + hInlet;
      end
      function [Dm, hOutlet] = process(comp,frequency,dInlet,pInlet,pOutlet,hInlet)
         Dm = flow(comp,frequency,dInlet);
         hOutlet = enthalpy(comp,pInlet,pOutlet,hInlet);
      end
   end
end