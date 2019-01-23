classdef Compressor   < matlab.mixin.Copyable
   properties
      Volume
      IsentropicEfficiency
      Dm
      DDm
      record
      h
      Tau
      pInlet
      dInlet
      pOutlet
      hInlet
   end
   methods
      function flow(comp,frequency,dInlet)
         Dm = 2*pi*50*frequency*dInlet*comp.Volume;
         comp.DDm = (Dm - comp.Dm)/comp.Tau;
      end
      function enthalpy(comp,pInlet,pOutlet,hInlet)
         s = CoolProp.PropsSI('S','H',hInlet,'P',pInlet,'CO2');
         hIdeal = CoolProp.PropsSI('H','S',s,'P',pOutlet,'CO2');
         comp.h = (hIdeal - hInlet)/comp.IsentropicEfficiency + hInlet;
      end
      function DDm = process(comp,t,x,frequency)
          % Inputs
%           dInlet = Inputs.dInlet;
%           pInlet = Inputs.pInlet;
%           pOutlet = Inputs.pOutlet;
%           hInlet = Inputs.hInlet;
%           frequency = Inputs.frequency;
          % State
          comp.Dm = x;
          % Process
          comp.flow(frequency,comp.dInlet);
          comp.enthalpy(comp.pInlet,comp.pOutlet,comp.hInlet);
          % Derivatives
          DDm = comp.DDm;
      end
      function initialize(comp,IsentropicEfficiency,Volume,Tau,Initial)
          comp.IsentropicEfficiency = IsentropicEfficiency;
          comp.Volume = Volume;
          comp.Tau = Tau;
          comp.Dm = Initial;
          % Record
          Record.t = [];
          Record.x = [];
          comp.record = Record;
      end
      function reinitialize(comp,Initial)
          comp.Dm = Initial;
          comp.record.t = [];
          comp.record.x = [];
      end
   end
end