classdef Comp
   properties
       Name
       E1 
       E2 
       G12  
       V12  
       tply  
       V21
       S
       Q
       zbar
       plystrain
       plystress
   end
   methods
       function obj = Comp(Name,E1,E2,G12,V12,tply)
          obj.Name = Name;
          obj.E1 = E1;
          obj.E2 = E2;
          obj.G12 = G12;
          obj.V12 = V12;
          obj.tply = tply;
          obj.V21= (V12*(E2/E1));
          obj.S = [1/obj.E1 -obj.V12*(obj.E2/obj.E1)/obj.E2 0;-obj.V12/obj.E1 1/obj.E2 0;0 0 1/obj.G12];
          obj.Q = [obj.E1/(1-obj.V12*obj.V21)   obj.V12*obj.E2/(1-obj.V12*obj.V21) 0;
                   obj.V12*obj.E2/(1-obj.V12*obj.V21)   obj.E2/(1-obj.V12*obj.V21) 0;
                   0    0   obj.G12];
       end
       function obj = rotate(obj,theta)
           m=cosd(theta);
           n=sind(theta);
           Ts=[m^2 n^2 2*m*n; n^2 m^2 -2*m*n; -m*n m*n m^2-n^2];
           Te=[m^2 n^2 m*n; n^2 m^2 -m*n; -2*m*n 2*m*n m^2-n^2];
           obj.Q=Ts\obj.Q*Te;
       end
   end
end