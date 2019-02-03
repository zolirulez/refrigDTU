classdef Fan < matlab.mixin.Copyable
   properties
       DV
       DDV
       d
       T
       maxVolumeFlow
       Tau
       record
   end
   methods
      function flow(fan,capacityRatio)
         DV = capacityRatio*fan.maxVolumeFlow;
         fan.DDV = (DV - fan.DV)/fan.tau;
      end
      function DDm = process(fan,t,x,capacityRatio)
          % State
          fan.Dm = x;
          % Process
          fan.flow(capacityRatio);
          % Derivatives
          DDm = fan.DDm;
      end
      function initialize(fan,maxVolumeFlow,Tau,Initial)
          fan.maxVolumeFlow = maxVolumeFlow;
          fan.Tau = Tau;
          fan.DV = Initial.DV;
          fan.T = Initial.T;
          % Record
          Record.t = [];
          Record.x = [];
          fan.record = Record;
      end
      function reinitialize(fan,Initial)
          fan.Dm = Initial.Dm;
          fan.T = Initial.T;
          fan.record.t = [];
          fan.record.x = [];
      end
      function pT2d(fan)
          fan.d = CoolProp.PropsSI('D','P',1.013e5,'T',fan.T,'AIR');
      end
   end
end