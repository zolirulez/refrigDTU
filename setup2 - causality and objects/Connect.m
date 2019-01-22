classdef Connect   < matlab.mixin.Copyable
   properties
   end
   methods
      function inlet(con,Output,Input,States)
          for it = 1:length(States)
              eval(['Output.' States{it} 'Inlet ='...
                  'Input.' States{it} '(end);']);
          end
      end
      function outlet(con,Output,Input,States)
          for it = 1:length(States)
              eval(['Output.' States{it} 'Outlet ='...
                  'Input.' States{it} '(1);']);
          end
      end
   end
end