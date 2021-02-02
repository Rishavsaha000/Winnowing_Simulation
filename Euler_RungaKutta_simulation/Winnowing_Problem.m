%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Winnowing Problem script%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                 %%%%%%%%%%Author - RISHAV SAHA - 229604 %%%%%%%%%%%%%%%
% Solving of a system of differential equations with Euler or forth order Runge-Kutta    
clc;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Solver selection%%%%%%%%%%%%%%%%%%%%%
Solver = input("Select a Solver: Press 1 for Euler or press 2 for RK4 -") ; %1 = Euler and 2 = RK4

%%%%%%%%%%%%%%%%%%%Fluid properties%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rho_F = 1.2;                           %kg/m3
Mu_F = 1.8*10^-5;                      %Pa.s

%%%%%%%%%%%%%%%%%Heavy particles 1 (Grains)Properties %%%%%%%%%%
Rho_P1 = 750;                           %kg/m3
Dp1 = 2.5*(10^-3);                      %m
Vp1 = (pi*Dp1^3)/6;                     %m3

%%%%%%%%%%%%%%%%%%Light Particles 2 (Chaff) Properties%%%%%%%%%%
Rho_P2 = 50;                           %kg/m3
Dp2 = 3.25*(10^-3);                    %m
Vp2 = (pi*Dp2^3)/6;                    %m3

%%%%%%%%%%%%%%%%Known values%%%%%%%%%%%%%%%%%%
g = 9.81;                              %m/s2, accleration due to gravity
h = 0.1;               	               %m, Height of rectangular slot
u0 = 0.2;                              %m/s, jet velocity

% Start of the computation interval
tstart = 0;
% End of the computation interval
tend = 2;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Forces acting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_Buoyancy1 = -Rho_F*Vp1*g;                       %only on Y direction(Particle 1)
F_gravity1 = Rho_P1*Vp1*g;                        %only on Y direction(Particle 1)
F_Buoyancy2 = -Rho_F*Vp2*g;                        %only on Y direction(Particle 2)
F_gravity2 = Rho_P2*Vp2*g;                         %only on Y direction(Particle 2)

Rho_P1XVp1= Rho_P1*Vp1;
Rho_P2XVp2= Rho_P2*Vp2;




%%%%%%%%%%%%%%%%%%%%%%%%eulerian velocity field%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uFx1= @(x1,y1) 6.2*u0*((h/x1)^0.5)*exp(-50*((y1/x1)^2));          %For particle 1 
uFy1 = 0;
uFx2 = @(x2,y2) 6.2*u0*((h/x2)^0.5)*exp(-50*(y2/x2)^2);           %For particle 2
uFy2 = 0;

%inital position of the particles%
x0 = 0.5;
y0 = 0.5;
xc = 0.55;
yc = -0.5;
x3 = [xc xc];                                  %BinSeparator
y3 = [yc yc+0.05];                             %BinSeparator
x4 = [0 2];                                    %Bin
y4 = [yc yc];                                  %Bin
  
            %%%%%%%%%%%%%%%%%Solving ODEs both for both X and Y direction%%%%%%%%%%
            
            
%RHS of ODE equation 3.52 from the script for X and Y Axis
f1_tu_x = @(t,F_drag1x)  F_drag1x/Rho_P1XVp1;                           %for Particle 1 in X direction 
f1_tu_y = @(t,F_drag1y) (F_Buoyancy1+F_gravity1+ F_drag1y)/Rho_P1XVp1;  %for Particle 1 in Y direction

f2_tu_x = @(t,F_drag2x)  F_drag2x/Rho_P2XVp2;                           %for Particle 2 in X direction
f2_tu_y = @(t,F_drag2y) (F_Buoyancy2+F_gravity2+ F_drag2y)/Rho_P2XVp2;  %for Particle 2 in Y direction



%RHS of ODE equation 3.51 from the script for X and Y Axis
f1_tx_x = @(t,uP1x) uP1x;       %for Particle 1 in X direction
f1_ty_y = @(t,uP1y) uP1y;       %for Particle 1 in Y direction

f2_tx_x = @(t,uP2x) uP2x;       %for Particle 2 in X direction
f2_ty_y = @(t,uP2y) uP2y;       %for Particle 2 in Y direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Computation loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Number of time step differences for which the accuracy is calculated%%%%%
n = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%time arrays%%%%%%%%%%%%%%%%%%%%%%
dt = zeros(1,n);


if Solver == 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%EULER METHOD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for j = 1:n
    %calculate time step
    dt(j) = 10^(-j);
    
    % Allocate t
    t = tstart:dt(j):tend;

    %%%%%%%%%%%%%%%%%%%%%%Allocate the Solution vectors with initial point%%%%%%%%%%%%%%%%%%%%%%%%%
    x1 = ones(1,length(t))*x0;                %Intitially Particle 1 at 0.5m in X direction
    y1 = ones(1,length(t))*y0;                %Intitially Particle 1 at 0.5m in Y direction
    x2 = ones(1,length(t))*x0;                       %Intitially Particle 2 at 0.5m in X direction
    y2 = ones(1,length(t))*y0;                       %Intitially Particle 2 at 0.5m in Y direction

    uP1x = zeros(1,length(t));                 %Particle 1 velocity on X direction
    uP1y = zeros(1,length(t));                 %Particle 1 velocity on Y direction
    uP2x = zeros(1,length(t));                        %Particle 2 velocity on X direction
    uP2y = zeros(1,length(t));                        %Particle 2 velocity on Y direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Forces acting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_drag1x = zeros(1,length(t));                    %Particle 1 Drag force on X direction
    F_drag1y = zeros(1,length(t));                    %Particle 1 Drag force on Y direction
    F_drag2x = zeros(1,length(t));                     %Particle 2 Drag force on X direction
    F_drag2y = zeros(1,length(t));                     %Particle 2 Drag force on Y direction
    
    
        for i = 1:length(t)-1
            
           %  uP1 = (uP1x(i)^2+ uP1y(i)^2).^0.5; 
           %  uP2 = (uP2x(i)^2+ uP2y(i)^2).^0.5;   
             Abs_Urel1 = ((uFx1(x1(i),y1(i))-uP1x(i))^2 + (0-uP1y(i))^2)^0.5;
             Abs_Urel2 = ((uFx2(x2(i),y2(i))-uP2x(i))^2 + (0-uP2y(i))^2)^0.5;
             ReP1 = (Rho_F*Abs_Urel1*Dp1)./Mu_F;
             ReP2 = (Rho_F*Abs_Urel2*Dp2)./Mu_F;
             
             
             if ReP1<800
             CD1 = (24*(1+0.15*ReP1^0.687))/ReP1;
             else 
             CD1 = 0.44;
             end
             
             if ReP2<800
             CD2 = (24*(1+0.15*ReP2^0.687))/ReP2;
             else 
             CD2 = 0.44;
             end
     
             F_drag1x(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(uFx1(x1(i),y1(i))-uP1x(i))*Abs_Urel1;
             F_drag1y(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(0-uP1y(i))*Abs_Urel1;
             
             
             F_drag2x(i) = 0.5*pi*Dp2^2*Rho_F*CD2*(uFx2(x2(i),y2(i))-uP2x(i))*Abs_Urel2;
             F_drag2y(i) = 0.5*pi*Dp2^2*Rho_F*CD2*(0-uP2y(i))*Abs_Urel2;
   
             %solving Equation 3.52
             uP1x(i+1) = uP1x(i) + dt(j)*f1_tu_x(t(i),F_drag1x(i));
             uP1y(i+1) = uP1y(i) + dt(j)*f1_tu_y(t(i),F_drag1y(i)); 
             
             uP2x(i+1) = uP2x(i) + dt(j)*f2_tu_x(t(i),F_drag2x(i));
             uP2y(i+1) = uP2y(i) + dt(j)*f2_tu_y(t(i),F_drag2y(i));
    
             %solving Equation 3.51
             x1(i+1) = x1(i) + dt(j)*f1_tx_x(t(i),uP1x(i+1));
             y1(i+1) = y1(i) - dt(j)*f1_ty_y(t(i),uP1y(i+1));
             
             x2(i+1) = x2(i) + dt(j)*f2_tx_x(t(i),uP2x(i+1));
             y2(i+1) = y2(i) - dt(j)*f2_ty_y(t(i),uP2y(i+1));
        
                                          
        end
        
        if     j == 1
             figure("Name","Accuracy check of Euler's method");
             plot(x1,y1,x2,y2,"k:","linewidth",2);
        elseif j == 2
             plot(x1,y1,x2,y2,"b:","linewidth",2);
        elseif j == 3
             plot(x1,y1,x2,y2,"r:","linewidth",2);
        elseif j == 4
             plot(x1,y1,x2,y2,"y:","linewidth",2);
        elseif j == 5
             plot(x1,y1,x2,y2,"m:","linewidth",2);
             ylabel("Y- axis (Height [m])"); 
             xlabel("X- axis (Distance [m])");
             legend("for grains (dt=10^-1) ","for chaffs (dt=10^-1)"...
                 ,"for grains (dt=10^-2) ","for chaffs (dt=10^-2)"...
                 ,"for grains (dt=10^-3) ","for chaffs (dt=10^-3)"...
                 ,"for grains (dt=10^-4) ","for chaffs (dt=10^-4)"...
                 ,"most accurate for grains","most accurate for chaffs");
             axis([x0-0.05 x0+0.25 -0.5 y0]);
        end
        hold on
end
hold off
figure("Name","Trajectory using Euler's method");

elseif Solver ==2  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%fourth-order Runge-Kutta method%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
for j = 1:n
    %calculate time step
    dt(j) = 10^(-j);
    
    % Allocate t
    t = tstart:dt(j):tend;

    %%%%%%%%%%%%%%%%%%%%%%Allocate the Solution vectors with initial point%%%%%%%%%%%%%%%%%%%%%%%%%
    x1 = ones(1,length(t))*x0;                %Intitially Particle 1 at 0.5m in X direction
    y1 = ones(1,length(t))*y0;                %Intitially Particle 1 at 0.5m in Y direction
    x2 = ones(1,length(t))*x0;                       %Intitially Particle 2 at 0.5m in X direction
    y2 = ones(1,length(t))*y0;                       %Intitially Particle 2 at 0.5m in Y direction

    uP1x = zeros(1,length(t));                 %Particle 1 velocity on X direction
    uP1y = zeros(1,length(t));                 %Particle 1 velocity on Y direction
    uP2x = zeros(1,length(t));                        %Particle 2 velocity on X direction
    uP2y = zeros(1,length(t));                        %Particle 2 velocity on Y direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Forces acting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F_drag1x = zeros(1,length(t));                    %Particle 1 Drag force on X direction
    F_drag1y = zeros(1,length(t));                    %Particle 1 Drag force on Y direction
    F_drag2x = zeros(1,length(t));                     %Particle 2 Drag force on X direction
    F_drag2y = zeros(1,length(t));                     %Particle 2 Drag force on Y direction
  
             for i = 1:length(t)-1
                  uP1 = (uP1x(i)^2+ uP1y(i)^2).^0.5; 
                  uP2 = (uP2x(i)^2+ uP2y(i)^2).^0.5;   
                  Abs_Urel1 = ((uFx1(x1(i),y1(i))-uP1x(i))^2 + (0-uP1y(i))^2)^0.5;
                  Abs_Urel2 = ((uFx2(x2(i),y2(i))-uP2x(i))^2 + (0-uP2y(i))^2)^0.5;
                  ReP1 = (Rho_F*Abs_Urel1*Dp1)./Mu_F;
                  ReP2 = (Rho_F*Abs_Urel2*Dp2)./Mu_F;
             
                  
                  if ReP1<800
                  CD1 = (24*(1+0.15*ReP1^0.687))/ReP1;
                  else 
                  CD1 = 0.44;
                  end
                  
                  if ReP2<800
                  CD2 = (24*(1+0.15*ReP2^0.687))/ReP2;
                  else 
                  CD2 = 0.44;
                  end   
                  F_drag1x(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(uFx1(x1(i),y1(i))-uP1x(i))*Abs_Urel1;
                  F_drag1y(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(0-uP1y(i))*Abs_Urel1;  
             
                  F_drag2x(i) = 0.5*pi*Dp2^2*Rho_F*CD2*(uFx2(x2(i),y2(i))-uP2x(i))*Abs_Urel2;
                  F_drag2y(i) = 0.5*pi*Dp2^2*Rho_F*CD2*(0-uP2y(i))*Abs_Urel2;
    
                 %solving Equation 3.52   
    
                  k1x = f1_tu_x(t(i)             ,F_drag1x(i));
                  k2x = f1_tu_x(t(i) + dt(j)/2   ,F_drag1x(i)+(dt(j)/2)*k1x);
                  k3x = f1_tu_x(t(i) + dt(j)/2   ,F_drag1x(i)+(dt(j)/2)*k2x);
                  k4x = f1_tu_x(t(i) + dt(j)     ,F_drag1x(i)+ dt(j)*k3x);
                  uP1x(i+1) = uP1x(i) + (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x );
    
                  k1y = f1_tu_y(t(i) ,  F_drag1y(i));
                  k2y = f1_tu_y(t(i) + dt(j)/2,  F_drag1y(i)+   0.5*dt(j)*k1y);
                  k3y = f1_tu_y(t(i) + dt(j)/2,  F_drag1y(i)+   0.5*dt(j)*k2y);
                  k4y = f1_tu_y(t(i) + dt(j)  ,  F_drag1y(i) +  dt(j)*k3y);
                  uP1y(i+1) = uP1y(i) + (1/6)*dt(j)*( k1y + 2*k2y + 2*k3y + k4y  );
    
                  k1x = f2_tu_x(t(i),  F_drag2x(i));
                  k2x = f2_tu_x(t(i) + dt(j)/2,  F_drag2x(i)+(dt(j)/2)*k1x);
                  k3x = f2_tu_x(t(i) + dt(j)/2,  F_drag2x(i)+(dt(j)/2)*k2x);
                  k4x = f2_tu_x(t(i) + dt(j)  ,  F_drag2x(i)  + dt(j)*k3x);
                  uP2x(i+1) = uP2x(i) + (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x );
    
    
                  k1x = f2_tu_y(t(i) ,  F_drag2y(i));
                  k2x = f2_tu_y(t(i) + dt(j)/2,  F_drag2y(i)+(dt(j)/2)*k1x);
                  k3x = f2_tu_y(t(i) + dt(j)/2,  F_drag2y(i)+(dt(j)/2)*k2x);
                  k4x = f2_tu_y(t(i) + dt(j)  ,  F_drag2y(i)  + dt(j)*k3x);
                  uP2y(i+1) = uP2y(i) + (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x );              
                  
                  
                  %solving Equation 3.51
                  k1x = f1_tx_x(t(i)       ,  uP1x(i+1));
                  k2x = f1_tx_x(t(i) + dt(j)/2,  uP1x(i+1)+(dt(j)/2)*k1x);
                  k3x = f1_tx_x(t(i) + dt(j)/2,  uP1x(i+1)+(dt(j)/2)*k2x);
                  k4x = f1_tx_x(t(i) + dt(j)  ,  uP1x(i+1)+ dt(j)*k3x);
                  x1(i+1) = x1(i) + (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x  );
    
                  k1y = f1_ty_y(t(i),uP1y(i+1));
                  k2y = f1_ty_y(t(i) + dt(j)/2,  uP1y(i+1)+(dt(j)/2)*k1y);
                  k3y = f1_ty_y(t(i) + dt(j)/2,  uP1y(i+1)+(dt(j)/2)*k2y);
                  k4y = f1_ty_y(t(i) + dt(j)  ,  uP1y(i+1)+ dt(j)*k3y);
                  y1(i+1) = y1(i) - (1/6)*dt(j)*( k1y + 2*k2y + 2*k3y + k4y );
    

                  k1x = f2_tx_x(t(i)       ,  uP2x(i+1));
                  k2x = f2_tx_x(t(i) + dt(j)/2,  uP2x(i+1)+(dt(j)/2)*k1x);
                  k3x = f2_tx_x(t(i) + dt(j)/2,  uP2x(i+1)+(dt(j)/2)*k2x);
                  k4x = f2_tx_x(t(i) + dt(j)  ,  uP2x(i+1)+ dt(j)*k3x);
                  x2(i+1) = x2(i) + (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x );
    
                  k1x = f2_ty_y(t(i)       ,  uP2y(i+1));
                  k2x = f2_ty_y(t(i) + dt(j)/2,  uP2y(i+1)+(dt(j)/2)*k1x);
                  k3x = f2_ty_y(t(i) + dt(j)/2,  uP2y(i+1)+(dt(j)/2)*k2x);
                  k4x = f2_ty_y(t(i) + dt(j)  ,  uP2y(i+1)+ dt(j)*k3x);
                  y2(i+1) = y2(i) - (1/6)*dt(j)*( k1x + 2*k2x + 2*k3x + k4x );
             end  
        if    j == 1
             figure("Name","Accuracy check of RK4 method");
             plot(x1,y1,x2,y2,"k:","linewidth",2);
        elseif j == 2
             plot(x1,y1,x2,y2,"b:","linewidth",2);
        elseif j == 3
             plot(x1,y1,x2,y2,"r:","linewidth",2);
        elseif j == 4
             plot(x1,y1,x2,y2,"y:","linewidth",2);
        elseif j == 5
             plot(x1,y1,x2,y2,"m:","linewidth",2);
             ylabel("Y- axis (Height [m])"); 
             xlabel("X- axis (Distance [m])");
             legend("for grains (dt=10^-1) ","for chaffs (dt=10^-1)"...
                 ,"for grains (dt=10^-2) ","for chaffs (dt=10^-2)"...
                 ,"for grains (dt=10^-3) ","for chaffs (dt=10^-3)"...
                 ,"for grains (dt=10^-4) ","for chaffs (dt=10^-4)"...
                 ,"most accurate for grains","most accurate for chaffs");
             axis([x0-0.05 x0+0.25 -0.5 y0]);
        end
        hold on

end
figure("Name","Trajectory using RK4 method");

else
     disp('Please select a solver.')
     return
end

  %%%% Plots %%%
    
    grid on
    hold on 
    plot (x1,y1,'r:.',"linewidth",2);              %Plot for particle 1
    plot (x2,y2,'b:',"linewidth",2);               %Plot for particle 2
    plot (x3,y3,'k',"linewidth",3);                %Plot for BinSeparator
    plot (x4,y4,'k',"linewidth",3);
    plot (x0,y0,"*","linewidth",3);              %plot initial point for both Particles
    text(x0+0.005,y0,"Initial point");
    text(0.6,-0.45,"Bin 2");
    text(0.47,-0.45,"Bin 1");
    text(0.56,-0.02,"Air is blown horizontally");
    annotation('arrow', [0.6 0.7], [0.5 0.5]);
    ylabel("Y- axis (Height [m])"); 
    xlabel("X- axis (Distance [m])");
    legend("trajectory of grains","trajectory of chaffs","Bin");
    axis([x0-0.05 x0+0.15 -0.5 y0+0.2]);
    
    %%%%%%%%%END%%%%%%%%%%