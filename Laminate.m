classdef Laminate
   properties
      %% Defines the structure for Laminates 
      Name              %Name of the Laminate
      stackup           %Materials that make up the laminate
      orientation       %Direction of each material
      sym               %Is the material symmetric or not
      thickness = 0;    %Overall Thickness of the laminate
      ABD               %ABD Matrix
      A = zeros(3,3);                %A matrix
      B = zeros(3,3);                %B matrix
      D = zeros(3,3);                %D matrix
      abd               %abd matrix
   end
   methods 
       function obj = Laminate(Name,stackup,orientation,sym) % Function that defines a laminate
        %Size Check       
        if length(stackup)~=length(orientation) % Checks the size of the s
           disp('ERROR: Stackup and orientation matrix do not match')
           return
        end
       %Given from inputs
       obj.Name = Name;
       obj.stackup = stackup;           %Gives the object the stackup matrix
       obj.orientation = orientation;   %Gives the object the orientation matrix
       obj.sym = sym;                   %If symmetric the composite is automatically doubled and mirrored
       %Mirroring Process for Sym 
       if sym ==1
           layers = length(obj.stackup);
           for i = layers:-1:1
               obj.stackup(end+1)=obj.stackup(i);
               obj.orientation(end+1)=obj.orientation(i);
           end
           clear temp
       end
       %Thickness Calculatation 
       layers = length(obj.stackup);     %Number of layers in stackup
       for i = 1:layers
           obj.thickness=obj.thickness+obj.stackup(i).tply;
       end
       %ABD Caclulations
           for i= 1:layers              %Goes through each layer
           obj.stackup(i)=rotate(obj.stackup(i),obj.orientation(i)); %rotates the Q matrix using rotate function
           obj.A = obj.A + obj.stackup(i).Q*obj.stackup(i).tply;
           %Zbar calc
           even=rem(layers,2);           
           if even==0
              if i <=layers/2
                 temp2=0;
                 for j = i:layers/2
                    temp2=temp2+obj.stackup(j).tply;
                 end
                 obj.stackup(i).zbar = temp2-obj.stackup(i).tply/2;
              else
                 temp2=0;
                 for j = (layers/2+1):i
                    temp2=temp2+obj.stackup(j).tply;
                 end
                 obj.stackup(i).zbar = -(temp2-obj.stackup(i).tply/2);
              end
           else
               mt=obj.stackup(((layers-1)/2)+1).tply;
               if i <=(layers-1)/2
                 temp2=0;
                 for j = i:layers/2
                    temp2=temp2+obj.stackup(j).tply;
                 end
                 obj.stackup(i).zbar = mt+temp2-obj.stackup(i).tply/2;
               elseif i >=((layers-1)/2)+2
                 temp2=0;
                 for j = ((layers-1)/2+2):i
                    temp2=temp2+obj.stackup(j).tply;
                 end
                 obj.stackup(i).zbar = -(mt+temp2-obj.stackup(i).tply/2);
               else
                 obj.stackup(i).zbar =0;
               end            
           end
           %Calculate B and D from Z bar
           obj.B = obj.B + obj.stackup(i).Q*(obj.stackup(i).tply*obj.stackup(i).zbar);
           obj.D = obj.D + obj.stackup(i).Q*((obj.stackup(i).tply^3/12)+obj.stackup(i).tply*obj.stackup(i).zbar^2);
           end    
           obj.ABD=[obj.A,obj.B;obj.B,obj.D];
           obj.abd=inv(obj.ABD);
       end
       function printprop(obj) %Function to print the properities of a laminate
           disp('----------------------------------------------------------')
           fprintf('Stackup Sequence for Laminate: %s\n',obj.Name)
           disp('----------------------------------------------------------')
           temp='Material Name';
           temp2='Orientation';
           fprintf('%s%15s\n',temp,temp2)
           layers=length(obj.stackup);
           for i = 1:layers
               fprintf('%s\t\t\t%4.4g Degrees\n',obj.stackup(i).Name,obj.orientation(i))
           end
           disp('----------------------------------------------------------')
           disp('ABD Matrix')           
           for i = 1:6
               for j = 1:6
                    fprintf('%8.3g  ',obj.ABD(i,j))
               end
               fprintf('\n')
           end
           disp('----------------------------------------------------------')
           disp('A Matrix')           
           for i = 1:3
               for j = 1:3
                    fprintf('%10.5g  ',obj.ABD(i,j))
               end
               fprintf('\n')
           end
           disp('----------------------------------------------------------')
           disp('B Matrix')           
           for i = 1:3
               for j = 4:6
                    fprintf('%10.5g  ',obj.ABD(i,j))
               end
               fprintf('\n')
           end
           disp('----------------------------------------------------------')
           disp('D Matrix')           
           for i = 4:6
               for j = 4:6
                    fprintf('%10.5g  ',obj.ABD(i,j))
               end
               fprintf('\n')
           end
           disp('----------------------------------------------------------')
       end
       function obj = strainNM(obj,N,M)
           NM = [N;M];
           e0k= obj.abd*NM;
           layers = length(obj.stackup);
           for i = 1:layers               
               obj.stackup(i).plystrain=e0k(1:3)+obj.stackup(i).zbar*e0k(4:6);
               obj.stackup(i).plystress=obj.stackup(i).Q*obj.stackup(i).plystrain;
           end
           
       end
   end
end