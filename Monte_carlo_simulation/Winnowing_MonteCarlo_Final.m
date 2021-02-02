%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Winnowing Problem 4.6 script%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                 %%%%%%%%%%Author - RISHAV SAHA - 229604 %%%%%%%%%%%%%%%                                
clc; clear all;

%%%%%%%%%%%%%%%%%%%Fluid properties%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rho_F = 1.2;                           %kg/m3
Mu_F = 1.8*10^-5;                      %Pa.s

%%%%%%%%%%%%%%%%%Heavy particles 1 (Grains)Properties %%%%%%%%%%
Rho_P1 = 750;                          %kg/m3
Dp1 = 2.5*(10^-3);                     %m
Vp1 = (pi*Dp1^3)/6;                    %m3

%%%%%%%%%%%%%%%%%%Light Particles 2 (ChaffS) Properties%%%%%%%%%%
Rho_P2 = 50;                           %kg/m3
Dp2    = 3.25*(10^-3);                  %m
Vp2    = (pi*Dp2^3)/6;                  %m3

%%%%%%%%%%%%%%%%Known values%%%%%%%%%%%%%%%%%%
g = 9.81;                              %m/s2, accleration due to gravity
h = 0.1;               	               %m, Height of rectangular slot
u0 = 0.2;                              %m/s, jet velocity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Forces acting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_Buoyancy1 = -Rho_F*Vp1*g;                         %only on Y direction(Particle 1)
F_gravity1 =   Rho_P1*Vp1*g;                        %only on Y direction(Particle 1)
F_Buoyancy2 = -Rho_F*Vp2*g;                         %only on Y direction(Particle 2)
F_gravity2 =   Rho_P2*Vp2*g;                        %only on Y direction(Particle 2)

Rho_P1XVp1= Rho_P1*Vp1;
Rho_P2XVp2= Rho_P2*Vp2;

%%%%%%%%%%%%%%%%%%%%%%%%eulerian velocity field%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uFx = @(x,y)   6.2*u0*((h/x)^0.5)*exp(-50*((y/x)^2));  %only on X direction 

uP1x(1) = 0;                        %Particle 1 velocity on X direction
uP1y(1) = 0;                        %Particle 1 velocity on Y direction
uP2x(1) = 0;                        %Particle 2 velocity on X direction
uP2y(1) = 0;                        %Particle 2 velocity on Y direction

%inital position of the particles%
x0 = 0.5; y0 = 0.5;

%BinSeparator
xc = 0.55; yc = -0.5;
x3 = [xc xc];   
y3 = [yc yc+0.05];    

%Bin                      
x4 = [0 2];                                    
y4 = [yc yc];                                 

%RHS of ODE equation 3.52
f_tu_x = @(F_d,RhoXVp)          F_d/RhoXVp;             %X direction 
f_tu_y = @(F_d,F_B,F_g,RhoXVp) (F_B +F_g+ F_d)/RhoXVp;  %Y direction

%RHS of ODE equation 3.51
f_tx_x = @(uPx) uPx;       %X direction
f_ty_y = @(uPy) uPy;       %Y direction
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%selection%%%%%%%%%%%%%%%%%%%%%
Task = input("1=Accuracy gain   2=Trajectory of 1000 particles--") ;

if Task == 1

    %%%%%%%%%%%%%%%%%%%%%%%%Considered Accurate Reference Solution%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%Averaging Combined method 50 times%%%%%%%%%%%%%%%
    
    for k = 1:50 %Reference (Almost Accurate) solution
    N = 100; %number of Random samples

    x1_AU = x0;             
    y1_AU = y0;             
    x2_AU = x0;         
    y2_AU = y0;
    dt = 10^-4;
    i = 1;
       while y1_AU > -0.5  %Grains
            %solving Equation 3.51 BY EULER
             x1_AU(i+1) = x1_AU(i) + dt*uP1x(i);
             y1_AU(i+1) = y1_AU(i) - dt*uP1y(i); 
             
            %solving Equation 3.52 by combined Monte Carlo and Euler
              Sum_RHS_x = 0;
              Sum_RHS_y = 0;
              for j = 1:N       %For N random time steps
                dt_rand = dt*rand();                   
                x1_rand = x1_AU(i) + uP1x(i)*dt_rand ;
                y1_rand = y1_AU(i) - uP1y(i)*dt_rand ;
                
                Abs_U1_rand = ((uFx(x1_rand,y1_rand)-uP1x(i))^2 + (0-uP1y(i))^2)^0.5; %Absolute Up - Uf
                ReP1_rand = (Rho_F*Abs_U1_rand*Dp1)./Mu_F; 
                if ReP1_rand<800
                CD1_rand = (24*(1+0.15*ReP1_rand^0.687))/ReP1_rand;
                else 
                CD1_rand = 0.44;
                end     
             
                F_drag1x_rand = 0.5*pi*(Dp1^2)*Rho_F*CD1_rand*(uFx(x1_rand,y1_rand)-uP1x(i))*Abs_U1_rand;
                F_drag1y_rand = 0.5*pi*(Dp1^2)*Rho_F*CD1_rand*(0-uP1y(i))*Abs_U1_rand;
                
                f1_tu_x_rand =   F_drag1x_rand/Rho_P1XVp1;                           %for X direction 
                f1_tu_y_rand =  (F_Buoyancy1+F_gravity1+ F_drag1y_rand)/Rho_P1XVp1;  %for Y direction
                Sum_RHS_x = Sum_RHS_x + f1_tu_x_rand ;
                Sum_RHS_y = Sum_RHS_y + f1_tu_y_rand ; 
              end
              
             uP1x(i+1) = uP1x(i) + dt*Sum_RHS_x/N;
             uP1y(i+1) = uP1y(i) + dt*Sum_RHS_y/N;
             
             X1_AU(i,k) = x1_AU(i) ;   Y1_AU(i,k) = y1_AU(i);  
             i = i +1;
       end  %Grains computation ends
 i = 1;
       while y2_AU > -0.5  %Chaffs
        
            %solving Equation 3.51 BY EULER
             x2_AU(i+1) = x2_AU(i) + dt*uP2x(i);
             y2_AU(i+1) = y2_AU(i) - dt*uP2y(i); 
             
            %solving Equation 3.52 by combined Monte Carlo and Euler
              Sum_RHS_x = 0;
              Sum_RHS_y = 0;
              
              for j = 1:N                %For N random time steps
                dt_rand = dt*rand;                     
                x2_rand = x2_AU(i) + uP2x(i)*dt_rand ;
                y2_rand = y2_AU(i) - uP2y(i)*dt_rand ;
                              
                Abs_U2_rand = ((uFx(x2_rand,y2_rand)-uP2x(i))^2 + (0-uP2y(i))^2)^0.5; %Absolute Up - Uf
                ReP2_rand = (Rho_F*Abs_U2_rand*Dp2)./Mu_F; 
                if ReP2_rand<800
                CD2_rand = (24*(1+0.15*ReP2_rand^0.687))/ReP2_rand;
                else 
                CD2_rand = 0.44;
                end     
             
                F_drag2x_rand = 0.5*pi*(Dp2^2)*Rho_F*CD2_rand*(uFx(x2_rand,y2_rand)-uP2x(i))*Abs_U2_rand;
                F_drag2y_rand = 0.5*pi*(Dp2^2)*Rho_F*CD2_rand*(0-uP2y(i))*Abs_U2_rand;
                f2_tu_x_rand =   F_drag2x_rand/Rho_P2XVp2;                           %for Particle 1 in X direction 
                f2_tu_y_rand =  (F_Buoyancy2+F_gravity2+ F_drag2y_rand)/Rho_P2XVp2;  %for Particle 1 in Y direction

                Sum_RHS_x = Sum_RHS_x + f2_tu_x_rand ;
                Sum_RHS_y = Sum_RHS_y + f2_tu_y_rand ;
              end
             uP2x(i+1) = uP2x(i) + dt*Sum_RHS_x/N;
             uP2y(i+1) = uP2y(i) + dt*Sum_RHS_y/N; 
             
             X2_AU(i,k) = x2_AU(i) ;   Y2_AU(i,k) = y2_AU(i); 
             i = i +1;        
       end   %Chaffs Accurate computation ends
       
    end   
    
    X1_Accurate = sum(X1_AU,2)/k;
    X2_Accurate = sum(X2_AU,2)/k;
    Y1_Accurate = sum(Y1_AU,2)/k;
    Y2_Accurate = sum(Y2_AU,2)/k;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%Solutions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1:3 
    dt(k) = 10^-(k+1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Combined MC and EU%%%%%%%%%%%%%%%%%%%%%%%
    N =  [10 100];  %Number of Random samples
    for z = 1:length(N)
    %Allocate the Solution vectors with initial point%
    x1_MC = x0;     x1_EU = x0;            
    y1_MC = y0;     y1_EU = y0;           
    x2_MC = x0;     x2_EU = x0;      
    y2_MC = y0;     y2_EU = y0;
    
 i = 1;
       while y1_MC > -0.5  %Grains
            %solving Equation 3.51 BY EULER
             x1_MC(i+1) = x1_MC(i) + dt(k)*uP1x(i);
             y1_MC(i+1) = y1_MC(i) - dt(k)*uP1y(i); 
             
            %solving Equation 3.52 by combined Monte Carlo and Euler
              Sum_RHS_x = 0;
              Sum_RHS_y = 0;
              for j = 1:N(z)       %For N random time steps
                dt_rand = dt(k)*rand();                   
                x1_rand = x1_MC(i) + uP1x(i)*dt_rand ;
                y1_rand = y1_MC(i) - uP1y(i)*dt_rand ;
                
                Abs_U1_rand = ((uFx(x1_rand,y1_rand)-uP1x(i))^2 + (0-uP1y(i))^2)^0.5; %Absolute Up - Uf
                ReP1_rand = (Rho_F*Abs_U1_rand*Dp1)./Mu_F; 
                if ReP1_rand<800
                CD1_rand = (24*(1+0.15*ReP1_rand^0.687))/ReP1_rand;
                else 
                CD1_rand = 0.44;
                end     
             
                F_drag1x_rand = 0.5*pi*(Dp1^2)*Rho_F*CD1_rand*(uFx(x1_rand,y1_rand)-uP1x(i))*Abs_U1_rand;
                F_drag1y_rand = 0.5*pi*(Dp1^2)*Rho_F*CD1_rand*(0-uP1y(i))*Abs_U1_rand;
                
                f1_tu_x_rand =   F_drag1x_rand/Rho_P1XVp1;                           %for X direction 
                f1_tu_y_rand =  (F_Buoyancy1+F_gravity1+ F_drag1y_rand)/Rho_P1XVp1;  %for Y direction
                Sum_RHS_x = Sum_RHS_x + f1_tu_x_rand ;
                Sum_RHS_y = Sum_RHS_y + f1_tu_y_rand ; 
              end
              
             uP1x(i+1) = uP1x(i) + dt(k)*Sum_RHS_x/N(z);
             uP1y(i+1) = uP1y(i) + dt(k)*Sum_RHS_y/N(z);
             
             X1_MC(i,k) = x1_MC(i) ;   Y1_MC(i,k) = y1_MC(i);  
             i = i +1;
       end  %Grains computation ends
 i = 1;
       while y2_MC > -0.5  %Chaffs
        
            %solving Equation 3.51 BY EULER
             x2_MC(i+1) = x2_MC(i) + dt(k)*uP2x(i);
             y2_MC(i+1) = y2_MC(i) - dt(k)*uP2y(i); 
             
            %solving Equation 3.52 by combined Monte Carlo and Euler
              Sum_RHS_x = 0;
              Sum_RHS_y = 0;
              
              for j = 1:N(z)                %For N random time steps
                dt_rand = dt(k)*rand;                     
                x2_rand = x2_MC(i) + uP2x(i)*dt_rand ;
                y2_rand = y2_MC(i) - uP2y(i)*dt_rand ;
                              
                Abs_U2_rand = ((uFx(x2_rand,y2_rand)-uP2x(i))^2 + (0-uP2y(i))^2)^0.5; %Absolute Up - Uf
                ReP2_rand = (Rho_F*Abs_U2_rand*Dp2)./Mu_F; 
                if ReP2_rand<800
                CD2_rand = (24*(1+0.15*ReP2_rand^0.687))/ReP2_rand;
                else 
                CD2_rand = 0.44;
                end     
             
                F_drag2x_rand = 0.5*pi*(Dp2^2)*Rho_F*CD2_rand*(uFx(x2_rand,y2_rand)-uP2x(i))*Abs_U2_rand;
                F_drag2y_rand = 0.5*pi*(Dp2^2)*Rho_F*CD2_rand*(0-uP2y(i))*Abs_U2_rand;
                f2_tu_x_rand =   F_drag2x_rand/Rho_P2XVp2;                           %for Particle 1 in X direction 
                f2_tu_y_rand =  (F_Buoyancy2+F_gravity2+ F_drag2y_rand)/Rho_P2XVp2;  %for Particle 1 in Y direction

                Sum_RHS_x = Sum_RHS_x + f2_tu_x_rand ;
                Sum_RHS_y = Sum_RHS_y + f2_tu_y_rand ;
              end
             uP2x(i+1) = uP2x(i) + dt(k)*Sum_RHS_x/N(z);
             uP2y(i+1) = uP2y(i) + dt(k)*Sum_RHS_y/N(z); 
             
             X2_MC(i,k) = x2_MC(i) ;   Y2_MC(i,k) = y2_MC(i); 
             i = i +1;        
       end   %Chaffs computation ends
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%MC Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if z ==1
                 MC_err1_10 (k) =  (abs(X1_MC(end,k)-X1_Accurate(end)))/X1_Accurate(end);
                 MC_err2_10 (k) =  (abs(X2_MC(end,k)-X2_Accurate(end)))/X2_Accurate(end);    
            else 
                 MC_err1_100 (k) =  (abs(X1_MC(end,k)-X1_Accurate(end)))/X1_Accurate(end);
                 MC_err2_100 (k) =  (abs(X2_MC(end,k)-X2_Accurate(end)))/X2_Accurate(end);
            end 
    end   

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Euler Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 i = 1;
       while y1_EU > -0.5 %Grains
             Abs_Urel1 = ((uFx(x1_EU(i),y1_EU(i))-uP1x(i))^2 + (0-uP1y(i))^2)^0.5;       
             ReP1 = (Rho_F*Abs_Urel1*Dp1)./Mu_F;
             if ReP1<800
             CD1 = (24*(1+0.15*ReP1^0.687))/ReP1;
             else 
             CD1 = 0.44;
             end
             F_drag1x(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(uFx(x1_EU(i),y1_EU(i))-uP1x(i))*Abs_Urel1;
             F_drag1y(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(0-uP1y(i))*Abs_Urel1;
             
             %solving Equation 3.52
             uP1x(i+1) = uP1x(i) + dt(k)*f_tu_x(F_drag1x(i),Rho_P1XVp1);
             uP1y(i+1) = uP1y(i) + dt(k)*f_tu_y(F_drag1y(i),F_Buoyancy1,F_gravity1,Rho_P1XVp1); 
             
             %solving Equation 3.51
             x1_EU(i+1) = x1_EU(i) + dt(k)*f_tx_x(uP1x(i+1));
             y1_EU(i+1) = y1_EU(i) - dt(k)*f_ty_y(uP1y(i+1)); 
             
             
             X1_EU(i,k) = x1_EU(i);   Y1_EU(i,k) =  y1_EU(i) ; 
             i = i+1;                              
       end    %Grains computation ends
 i = 1;
       while y2_EU > -0.5 %Chaffs
             Abs_Urel2 = ((uFx(x2_EU(i),y2_EU(i))-uP2x(i))^2 + (0-uP2y(i))^2)^0.5;       
             ReP2 = (Rho_F*Abs_Urel2*Dp2)./Mu_F;
             if ReP2<800
             CD2 = (24*(1+0.15*ReP2^0.687))/ReP2;
             else 
             CD2 = 0.44;
             end
             F_drag2x(i) = 0.5*pi*(Dp2^2)*Rho_F*CD2*(uFx(x2_EU(i),y2_EU(i))-uP2x(i))*Abs_Urel2;
             F_drag2y(i) = 0.5*pi*(Dp2^2)*Rho_F*CD2*(0-uP2y(i))*Abs_Urel2;
             
             %solving Equation 3.52
             uP2x(i+1) = uP2x(i) + dt(k)*f_tu_x(F_drag2x(i),Rho_P2XVp2);
             uP2y(i+1) = uP2y(i) + dt(k)*f_tu_y(F_drag2y(i),F_Buoyancy2,F_gravity2,Rho_P2XVp2); 
            
             %solving Equation 3.51
             x2_EU(i+1) = x2_EU(i) + dt(k)*f_tx_x(uP2x(i+1));
             y2_EU(i+1) = y2_EU(i) - dt(k)*f_ty_y(uP2y(i+1)); 
             
             X2_EU(i,k) = x2_EU(i);     Y2_EU(i,k) =  y2_EU(i) ;
             
             i = i+1;                              
       end    %Chaffs computation ends
       
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%EU Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   EU_err1 (k) =  (abs(X1_EU(end,k)-X1_Accurate(end)))/X1_Accurate(end);
   EU_err2 (k) =  (abs(X2_EU(end,k)-X2_Accurate(end)))/X2_Accurate(end);    

end

%%%% Plots %%%%

figure("Name","Trajectory  of Particles dt = 10^-4");  
    grid on
    hold on
    plot (X1_MC(:,end),Y1_MC(:,end),'b',"linewidth",3);                %Plot for particle 1_MC
    plot (X2_MC(:,end),Y2_MC(:,end),'r',"linewidth",3);                %Plot for particle 2_MC
    plot (X1_EU(:,end),Y1_EU(:,end),'m:',"linewidth",3);               %Plot for particle 1_EU
    plot (X2_EU(:,end),Y2_EU(:,end),'y:',"linewidth",3);               %Plot for particle 2_EU
    plot (x3,y3,'k',"linewidth",3);                                    %Plot for BinSeparator
    plot (x4,y4,'k',"linewidth",3);
    plot (x0,y0,"*","linewidth",3);                                    %plot initial point for both Particles
    text (x0+0.005,y0,"Initial point");
    text (0.6,-0.45,"Bin 2");
    text (0.47,-0.45,"Bin 1");
    text (0.56,-0.02,"Air is blown horizontally");
    annotation('arrow', [0.6 0.7], [0.5 0.5]);
    ylabel("Y- axis (Height [m])"); 
    xlabel("X- axis (Distance [m])");
    legend("Grains(Combined Monte-Carlo)","Chaffs(Combined Monte-Carlo)","Grains(Euler's Method)","Chaffs((Euler's Method))","Bin");
    axis([x0-0.05 x0+0.15 -0.5 y0+0.2]);
    hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%For Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure("Name","Euler vs Accurate solution");
    hold on
    grid on
    loglog(dt,EU_err1,"r:","Marker",'*',"linewidth",2)
    loglog(dt,EU_err2,"b:","Marker",'o',"linewidth",2)
    legend("For Grains","For Chaffs");
    ylabel("Error"); 
    xlabel("time step (dt)");
 figure("Name","Combined Monte carlo vs Accurate solution");
    hold on
    grid on
    loglog(dt,MC_err1_10,"r:","Marker",'*',"linewidth",2)
    loglog(dt,MC_err1_100,"r","Marker",'*',"linewidth",2)
    loglog(dt,MC_err2_10,"b:","Marker",'o',"linewidth",2)
    loglog(dt,MC_err2_100,"b","Marker",'o',"linewidth",2)
    legend("10 random samples (Grains)","100 random samples (Grains)","10  random samples (Chaffs)","100  random samples (Chaffs)");
    ylabel("Error"); 
    xlabel("time step (dt)");
    %%%%%%%%END%%%%%%%%%%  

else 

%Number of Random PARTICLES falling
N =  1000;

Grain_Bin2 = 0;
Chaffs_Bin1 = 0;

for k = 1:N
    
%%%%%%%%%%%%%%%%%Heavy particles 1 (Grains)Properties %%%%%%%%%%
Rho_P1 = 750;                          %kg/m3
Dp1 = abs((2.5 + randn(1))*10^-3);                     %m
Vp1 = (pi*Dp1^3)/6;                    %m3

%%%%%%%%%%%%%%%%%%Light Particles 2 (ChaffS) Properties%%%%%%%%%%
Rho_P2 = abs(50 + 20*randn());                           %kg/m3
Dp2    = abs((2+3*rand())*(10^-3));                  %m
Vp2    = (pi*Dp2^3)/6;                  %m3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Forces acting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F_Buoyancy1 = -Rho_F*Vp1*g;                         %only on Y direction(Particle 1)
F_gravity1 =   Rho_P1*Vp1*g;                        %only on Y direction(Particle 1)
F_Buoyancy2 = -Rho_F*Vp2*g;                         %only on Y direction(Particle 2)
F_gravity2 =   Rho_P2*Vp2*g;                        %only on Y direction(Particle 2)

Rho_P1XVp1= Rho_P1*Vp1;
Rho_P2XVp2= Rho_P2*Vp2;

    %Allocate the Solution vectors with initial point%
       x1_EU = x0;            
       y1_EU = y0;           
       x2_EU = x0;      
       y2_EU = y0;
    
    dt = 10^-2; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% EU%%%%%%%%%%%%%%%%%%%%%%%
  i = 1;
       while y1_EU > -0.5 %Grains
             Abs_Urel1 = ((uFx(x1_EU(i),y1_EU(i))-uP1x(i))^2 + (0-uP1y(i))^2)^0.5;       
             ReP1 = (Rho_F*Abs_Urel1*Dp1)./Mu_F;
             if ReP1<800
             CD1 = (24*(1+0.15*ReP1^0.687))/ReP1;
             else 
             CD1 = 0.44;
             end
             F_drag1x(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(uFx(x1_EU(i),y1_EU(i))-uP1x(i))*Abs_Urel1;
             F_drag1y(i) = 0.5*pi*(Dp1^2)*Rho_F*CD1*(0-uP1y(i))*Abs_Urel1;
             
             %solving Equation 3.52
             uP1x(i+1) = uP1x(i) + dt*f_tu_x(F_drag1x(i),Rho_P1XVp1);
             uP1y(i+1) = uP1y(i) + dt*f_tu_y(F_drag1y(i),F_Buoyancy1,F_gravity1,Rho_P1XVp1); 
             
             %solving Equation 3.51
             x1_EU(i+1) = x1_EU(i) + dt*f_tx_x(uP1x(i+1));
             y1_EU(i+1) = y1_EU(i) - dt*f_ty_y(uP1y(i+1)); 
             
             
             X1_EU(i,k) = x1_EU(i);   Y1_EU(i,k) =  y1_EU(i) ; 
             
             i = i+1;                              
       end    %Grains computation ends
       
       if x1_EU(i) > xc
       Grain_Bin2 = Grain_Bin2 +1;
       else
       end
       
 i = 1;
       while y2_EU > -0.5 %Chaffs
             Abs_Urel2 = ((uFx(x2_EU(i),y2_EU(i))-uP2x(i))^2 + (0-uP2y(i))^2)^0.5;       
             ReP2 = (Rho_F*Abs_Urel2*Dp2)./Mu_F;
             if ReP2<800
             CD2 = (24*(1+0.15*ReP2^0.687))/ReP2;
             else 
             CD2 = 0.44;
             end
             F_drag2x(i) = 0.5*pi*(Dp2^2)*Rho_F*CD2*(uFx(x2_EU(i),y2_EU(i))-uP2x(i))*Abs_Urel2;
             F_drag2y(i) = 0.5*pi*(Dp2^2)*Rho_F*CD2*(0-uP2y(i))*Abs_Urel2;
             
             %solving Equation 3.52
             uP2x(i+1) = uP2x(i) + dt*f_tu_x(F_drag2x(i),Rho_P2XVp2);
             uP2y(i+1) = uP2y(i) + dt*f_tu_y(F_drag2y(i),F_Buoyancy2,F_gravity2,Rho_P2XVp2); 
            
             %solving Equation 3.51
             x2_EU(i+1) = x2_EU(i) + dt*f_tx_x(uP2x(i+1));
             y2_EU(i+1) = y2_EU(i) - dt*f_ty_y(uP2y(i+1)); 
             
             X2_EU(i,k) = x2_EU(i);     Y2_EU(i,k) =  y2_EU(i) ;
             
             i = i+1;                              
       end    %Chaffs computation ends
       
       if x2_EU(i) < xc
       Chaffs_Bin1 = Chaffs_Bin1 +1;
       else
       end
end


X1_EU(X1_EU==0) = NaN;    
X2_EU(X2_EU==0) = NaN;   
Y1_EU(Y1_EU==0) = NaN;   
Y2_EU(Y2_EU==0) = NaN;
Y1_EU(Y1_EU>0.5) = 0;   
Y2_EU(Y2_EU>0.5) = 0;

Grain_Bin2

Chaffs_Bin1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%selection%%%%%%%%%%%%%%%%%%%%%
Plot = input("1=Don't show plots   2=show plots--") ;
         if Plot == 2
            %%%% Plots %%%%
            figure("Name","Trajectory  of Grain Particles(dt = 10^-2 & 1000 samples)");  
                plot (X1_EU(:,:),Y1_EU(:,:),'b',"linewidth",1); 
                hold on
                plot (x3,y3,'k',"linewidth",3);                                    %Plot for BinSeparator
                plot (x4,y4,'k',"linewidth",3);
                plot (x0,y0,"*","linewidth",3);                                    %plot initial point for both Particles
                text (x0+0.005,y0,"Initial point");
                text (0.6,-0.45,"Bin 2");
                text (0.47,-0.45,"Bin 1");
                text (0.56,-0.02,"Air is blown horizontally");
                annotation('arrow', [0.6 0.7], [0.5 0.5]);
                ylabel("Y- axis (Height [m])"); 
                xlabel("X- axis (Distance [m])");
                axis([x0-0.05 x0+0.15 -0.5 y0+0.2]);
                hold off
              figure("Name","Trajectory  of Chaffs Particles(dt = 10^-2 & 1000 samples)");  
                plot (X2_EU(:,:),Y2_EU(:,:),'y',"linewidth",1);
                hold on
                plot (x3,y3,'k',"linewidth",3);                                    %Plot for BinSeparator
                plot (x4,y4,'k',"linewidth",3);
                plot (x0,y0,"*","linewidth",3);                                    %plot initial point for both Particles
                text (x0+0.005,y0,"Initial point");
                text (0.6,-0.45,"Bin 2");
                text (0.47,-0.45,"Bin 1");
                text (0.56,-0.02,"Air is blown horizontally");
                annotation('arrow', [0.6 0.7], [0.5 0.5]);
                ylabel("Y- axis (Height [m])"); 
                xlabel("X- axis (Distance [m])");
                axis([x0-0.05 x0+0.15 -0.5 y0+0.2]);
                hold off
            else
         end
end  