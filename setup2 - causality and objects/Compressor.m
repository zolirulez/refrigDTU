classdef Compressor   < matlab.mixin.Copyable
   properties
      Volume
      IsentropicEfficiency
      DmState
      DDm
      record
      hOutlet
   end
   methods
      function flow(comp,frequency,dInlet)
         Dm = 2*pi*frequency*dInlet*comp.Volume;
         comp.DDm = (Dm - comp.DmState)/comp.tau;
      end
      function enthalpy(comp,pInlet,pOutlet,hInlet)
         s = CoolProp.PropsSI('S','H',hInlet,'P',pInlet,'CO2');
         hOutletIdeal = CoolProp.PropsSI('H','S',s,'P',pOutlet,'CO2');
         comp.hOutlet = (hOutletIdeal - hInlet)/comp.IsentropicEfficiency + hInlet;
      end
      function DDm = process(comp,t,x,Inputs)
          % Inputs
          dInlet = Inputs.dInlet;
          pInlet = Inputs.pInlet;
          pOutlet = Inputs.pOutlet;
          hInlet = Inputs.hInlet;
          frequency = Inputs.frequency;
          % State
          comp.DmState = x;
          % Process
          comp.flow(comp,frequency,dInlet);
          comp.enthalpy(comp,pInlet,pOutlet,hInlet);
          % Derivatives
          DDm = comp.DDm;
      end
      function initialize(comp,IsentropicEfficiency,Volume,Tau,Initial)
          comp.IsentropicEfficiency = IsentropicEfficiency;
          comp.Volume = Volume;
          comp.tau = Tau;
          comp.DmState = Initial;
          % Record
          Record.t = [];
          Record.x = [];
          comp.record = Record;
      end
      function reinitialize(comp,Initial)
          comp.DmState = Initial;
          comp.record.t = [];
          comp.record.x = [];
      end
   end
end