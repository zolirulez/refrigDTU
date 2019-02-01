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
      function flow(comp,frequency)
         Dm = 2*pi*50*frequency*comp.dInlet*comp.Volume;
         comp.DDm = (Dm - comp.Dm)/comp.Tau;
      end
      function enthalpy(comp)
         s = CoolProp.PropsSI('S','H',comp.hInlet,'P',comp.pInlet,'CO2');
         try
             hIdeal = CoolProp.PropsSI('H','S',s,'P',comp.pOutlet,'CO2');
         catch
             hIdeal = 525e3; % CLEAR THIS TODO
         end
         comp.h = (hIdeal - comp.hInlet)/comp.IsentropicEfficiency + comp.hInlet;
      end
      function DDm = process(comp,t,x,frequency)
          % State
          comp.Dm = x;
          % Process
          comp.flow(frequency);
          %comp.enthalpy();
          % Derivatives
          DDm = comp.DDm;
      end
      function initialize(comp,IsentropicEfficiency,Volume,Tau,Initial)
          comp.IsentropicEfficiency = IsentropicEfficiency;
          comp.Volume = Volume;
          comp.Tau = Tau;
          comp.Dm = Initial.Dm;
          comp.h = Initial.h;
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