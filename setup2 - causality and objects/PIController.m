classdef PIController < matlab.mixin.Copyable
   properties
       err  % Error
       int  % Integrated error
       mn   % Minimum output
       mx   % Maximum output
       Ts   % Sampling time
       P    % Proportional gain
       I    % Integral gain
       neg  % Invert feedback or not (-1 for inversion or 1)
   end
   methods
      function out = react(ctrl,ref,meas)
          ctrl.err = ctrl.neg*(ref - meas);
          ctrl.int = ctrl.err*ctrl.Ts + ctrl.int;
          out = min(max(ctrl.P*ctrl.err + ctrl.I*ctrl.int,ctrl.mn),ctrl.mx);
      end
      function initialize(ctrl,K,Ti,initialIntegration,minOutput,maxOutput,negation,Ts)
          ctrl.Ts = Ts;
          ctrl.P = K;
          ctrl.I = ctrl.P/Ti;
          ctrl.mn = minOutput;
          ctrl.mx = maxOutput;
          ctrl.neg = negation;
          ctrl.int = initialIntegration;
      end
   end
end