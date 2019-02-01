classdef Fan < matlab.mixin.Copyable
   properties
      Kv
      tau
      Dm
      DDm
      record
      h
      pInlet
      pOutlet
      hInlet
      dInlet
   end
   methods
      function flow(valve,capacityRatio)
         % Kv value is in non-SI units!
         deltap = max([0,valve.pInlet - valve.pOutlet]);
         Dm = 8.7841e-06*capacityRatio*valve.Kv*sqrt(valve.dInlet*(deltap)) ;
         valve.DDm = (Dm - valve.Dm)/valve.tau;
      end
      function enthalpy(valve)
          valve.h = valve.hInlet;
      end
      function DDm = process(valve,t,x,capacityRatio)
          % State
          valve.Dm = x;
          % Process
          valve.flow(capacityRatio);
          %valve.enthalpy();
          % Derivatives
          DDm = valve.DDm;
      end
      function initialize(valve,Kv,Tau,Initial)
          valve.Kv = Kv;
          valve.tau = Tau;
          valve.Dm = Initial.Dm;
          valve.h = Initial.h;
          % Record
          Record.t = [];
          Record.x = [];
          valve.record = Record;
      end
      function reinitialize(valve,Initial)
          valve.Dm = Initial;
          valve.record.t = [];
          valve.record.x = [];
      end
   end
end