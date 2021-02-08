%% ELEC 4700- Assignment 1
% Saifuddin Mohammed, #101092039

% P1- Electron Modelling


set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth',2); 
close all

%Given Paramters
T = 300;  %Semicondctor Temperature 
Cmo = 9.10938356e-31; %Rest mass of the Electron
Cm = 0.26*Cmo; %Given Effective Mass of Electrons


% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt

Vth = sqrt(2*K*T/Cm) %Thermal Velocity Equation 


%Calculation of Mean Free Path, d
d = Vth*0.2e-12
d; 


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
x_R = 200e-9;
y_R = 100e-9;

%Setting Step size
st_size = 1e-9; 


electron_Population = 10000; %Assigned Population Size is 1000
%electron_Population = 10000; %Max

electron_Num = 50; %Electron population to be plotted

Tstep = y_R/Vth/100;
iter = 1000;
animation_plot = 0;


%Assigning all the Electrons positions

post = zeros(electron_Population, 4);

Temp = zeros(iter,1);


Trj = zeros(iter, electron_Num*2);






for ni = 1:electron_Population     %Producing a population of electrons.
    
    theta = rand*2*pi;            %Positional angle determination 
    
    post(ni,:) = [x_R*rand y_R*rand Vth*cos(theta) Vth*sin(theta)];   %The location of the electrons
    
    
end


%Utilizng time step and given paramters to make sure particles bounce off
%Storing and updating the position after every iteration

for ni = 1:iter
    post(:,1:2) = Tstep.*post(:,3:4) + post(:,1:2) ;

    
    %Boundary Conditions
    nk = post(:,1) > x_R;               %Collision detection 
    
    post(nk,1) = post(nk,1) - x_R;

    nk = post(:,1) < 0;
    post(nk,1) = post(nk,1) + x_R;

    nk = post(:,2) > y_R;
    post(nk,2) = 2*y_R - post(nk,2);
    post(nk,4) = -post(nk,4);

    nk = post(:,2) < 0;
    post(nk,2) = -post(nk,2);
    post(nk,4) = -post(nk,4);

    Temp(ni) = (sum(post(:,3).^2) + sum(post(:,4).^2))*Cm/K/2/electron_Population;

    
    %Keep track of trajectories
    
    for nk=1:electron_Num
        
        Trj(ni, (nk*2):(1+2*nk)) =  post(nk, 1:2);
        
    end 

    
    if animation_plot && mod(ni,50) == 0
        figure(1);
        plot(2,1,1);
        
        hold off;
        plot(post(1:electron_Num,1)./st_size, post(1:electron_Num,2)./st_size, 'o');
        set(gca,'Color', [1 1 1]);
        axis([0 x_R/st_size 0 y_R/st_size]);
        title(sprintf('Electron Trajectories',electron_Num, electron_Population));
        xlabel('Xp (nm)');
        ylabel('Yp (nm)');
        
        if ni > 1        %Semiconductor Temperature Plot Generation 
            figure(2); 
            plot(2,1,1);
            hold off;
            plot(Tstep*(0:ni-1), Temp(1:ni));
            set(gca,'Color', [0 0 0]);
            a_x = gca; 
           
            a_x.GridAlpha = 0.5;  % Make grid lines less transparent.
            a_x.GridColor = [1, 1, 1]; % Dark Green.

            axis([0 Tstep*iter min(Temp)*0.98 max(Temp)*1.02]);
            
            title('Temperature Plot');
            xlabel('T (s)');
            ylabel('Temperature (K)');
            grid on; 
        end
        
        pause(0.05);
        
    end
end


% Electron Trajectory represenation after the completion of Animation 

figure(1);
subplot(2,1,1);
title(sprintf('Electron Trajectories, SM1010192039',electron_Num, electron_Population));

xlabel(' Xp(nm)');

ylabel('Yp (nm)');

axis([0 x_R/st_size 0 y_R/st_size]);     %Setting axis to record values 

hold on;

for ni=1:electron_Num        %Storing the trajectory after completion 
    
    plot(Trj(:,ni*2)./st_size, Trj(:,ni*2+1)./st_size, '.');
    set(gca,'Color', [1 1 1]);
    
end

if(~animation_plot)
    
    figure(2);
    subplot(2,1,1);

    hold off;
    plot(Tstep*(0:iter-1), Temp);
    set(gca,'Color', [0 0 0]);
     a_x = gca; 
     a_x.GridAlpha = 0.2;  
     a_x.GridColor = [1, 1, 1]; 

    axis([0 Tstep*iter min(Temp)*0.92 max(Temp)*1.01]);
    
    
    title('Silicon Temperature, SM101092039');
    xlabel('T(s)');
    ylabel('Temperature(K)');
    grid on;
end



% P2 - Collisions with Mean Free Path (MFP) 


%Given Paramters
T = 300;  %Semicondctor Temperature 
Cmo = 9.10938356e-31; %Rest mass of the Electron
Cm = 0.26*Cmo; %Given Effective Mass of Electrons


% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt

Vth = sqrt(2*K*T/Cm) %Thermal Velocity Equation 


%Calculation of Mean Free Path, d
d = Vth*0.2e-12
d; 


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
x_R = 200e-9;
y_R = 100e-9;


%Setting the step size
st_size = 1e-9;

electron_Population = 10000; %Assigned Population Size is 1000
%electron_Population = 10000; %Max


%Using Maxewell-Boltzmann Distribution to Generate random velocities
Distr_MB = makedist('Normal', 0, sqrt(K*T/Cm));


%Initialization of the Position of Random Particles
%Setting up x and positions

for ni = 1:electron_Population
    
    theta = rand*2*pi;            %positional angle determination
    
    post(ni,:) = [x_R*rand y_R*rand random(Distr_MB) random(Distr_MB)];
end

%Calculation of the average Velocity of the distribution 

Vel_Avg = sqrt(sum(post(:,3).^2)/electron_Population + sum(post(:,4).^2)/electron_Population)

figure(3)                       %Histogram of Average velocities with Maxewell Boltzmann Distribution  
subplot(3,1,1);
z = sqrt(post(:,3).^2 + post(:,4).^2);

title('Electron Velocity Distribution Histogram');
histogram(z);
title('Electron Velocity Distribution Histogram, SM101092039'); 
xlabel('Velocity (m/s)');
ylabel('Total Electrons');


%Exponential Scattering Probability equation
P_scat = 1 - exp(-Tstep/0.2e-12)


for ni = 1:iter   % Performs the simulation with Random Velocities with updating positions. 
    
    
    %Boundary Conditions 
    post(:,1:2) = post(:,1:2) + Tstep.*post(:,3:4);

    nk = post(:,1) > x_R;
    post(nk,1) = - x_R + post(nk,1); 

    nk = post(:,1) < 0;
    post(nk,1) = x_R + post(nk,1); 

    nk = post(:,2) > y_R;
    post(nk,2) = 2*y_R - post(nk,2);
    post(nk,4) = -post(nk,4);

    nk = post(:,2) < 0;
    post(nk,2) = -post(nk,2);
    post(nk,4) = -post(nk,4);

   
    nk = rand(electron_Population, 1) < P_scat;
    post(nk,3:4) = random(Distr_MB, [sum(nk),2]);

   
    
    Temp(ni) = (sum(post(:,3).^2) + sum(post(:,4).^2))*Cm/K/2/electron_Population; % Tracking the Temperature of the Semicondctor 

    
  %Tracking the location of the particles for updating
  for nk=1:electron_Num
        
       Trj(ni, (2*nk):(2*nk+1)) = post(nk, 1:2);
        
  end 
    

    
    if animation_plot && mod(ni,10) == 0      %Tracking and stroing the data after iteration
       
        figure(4);                      %Plotting the Particle Trajectories                 
        plot(3,1,1);
        hold off;
        
        plot(post(1:electron_Num,1)./st_size, post(1:electron_Num,2)./st_size, 'o');
        axis([0 x_R/st_size 0 y_R/st_size]);
        title(sprintf('Electron Trajectories', electron_Num, electron_Population));
    
        xlabel('Xp (nm)');
        ylabel('Yp (nm)');
        
        
        if ni > 1              %Plotting the Temperature variation of the Silicon Semicondctor
            figure(5)
            plot(3,1,1);
            hold off;
            plot(Tstep*(0:ni-1), Temp(1:ni));
            
             set(gca,'Color', [0 0 0]);
            a_x = gca; 
           
            a_x.GridAlpha = 0.5;  
            a_x.GridColor = [1, 1, 1];
            
            axis([0 Tstep*iter min(Temp)*0.92 max(Temp)*1.01]);
            title('Silicon Temperature, SM101092039');
            xlabel('Time (s)');
            ylabel('Temperature (K)');
            grid on; 
        end

        
        
    end
end


%Final Plot after the completion of the animation

figure(4);                    %Plotting the final trajectories once the iterations have been completed
subplot(3,1,1);

title(sprintf('Electron Trajectories, SM101092039',electron_Num, electron_Population));
xlabel('Xp (nm)');
ylabel('Yp (nm)');
axis([0 x_R/st_size 0 y_R/st_size]);
hold on;


for ni=1:electron_Num           %Storing the Trajectory after completion 
    
    plot(Trj(:,ni*2)./st_size, Trj(:,ni*2+1)./st_size, '.');
    
    
end




if (~animation_plot)       %Final Temperature variance plot after complete iteration 
    
    figure(5)
    subplot(3,1,1);
    hold off;
    plot(Tstep*(0:iter-1), Temp);
    
   set(gca,'Color', [0 0 0]);
    a_x = gca; 
    a_x.GridAlpha = 0.5;  
    a_x.GridColor = [1, 1, 1];
    
    axis([0 Tstep*iter min(Temp)*0.92 max(Temp)*1.02]);
    title('Silicon Temperature, SM101092039');
    xlabel('T(s)');
    ylabel('Temperature (K)');
    grid on;
    
end


% P3 - Enhancements 



%Given Paramters
T = 300;  %Semicondctor Temperature 
Cmo = 9.10938356e-31; %Rest mass of the Electron
Cm = 0.26*Cmo; %Given Effective Mass of Electrons


% Calculation of Thermal Velocity
K = 1.38064852e-23;  %Boltzmann Constatnt

Vth = sqrt(2*K*T/Cm) %Thermal Velocity Equation 


%Calculation of Mean Free Path, d
d = Vth*0.2e-12
d; 


%Given Nominal Dimensions of Semiconductor 200 nm x 100 nm
x_R = 200e-9;
y_R = 100e-9;


%Setting the step size
st_size = 1e-9;



electron_Num = 10;

% If test_conditon = 0, then Diffusion
test_condition =1; 


% Determination of Rectangle Boxes
boxes = st_size.*[80 120 0 40; 80 120 60 100];
spec_boxes = [0 1];


%Initialization of the Position of Random Particles
%Setting up x and positions
for ni = 1:electron_Population
    
    theta = rand*2*pi;            %positional angle determination
    
    post(ni,:) = [x_R*rand y_R*rand random(Distr_MB) random(Distr_MB)];

   
    
    while(box_no(post(ni,1:2), boxes))  %Using the function box_no to prevent the present of electrons within the boxes
        post(ni,1:2) = [x_R*rand y_R*rand];
    end
end



%Performing simulation with Random velocitieis and updating the positons
%regularly 

for ni = 1:iter
    
    %Defining the boundary conditions
    
    post(:,1:2) = post(:,1:2) + Tstep.*post(:,3:4);   %New position calculation 

    nk = post(:,1) > x_R;           
    post(nk,1) = post(nk,1) - x_R;

    nk = post(:,1) < 0;
    post(nk,1) = post(nk,1) + x_R;

    nk = post(:,2) > y_R;

    if(test_condition)
        post(nk,2) = 2*y_R - post(nk,2);
        post(nk,4) = -post(nk,4);
        
    else
        % Diffusive Condition 
        
        post(nk,2) = y_R;
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = -abs(Z.*sin(theta));
        
    end

    
    nk = post(:,2) < 0;

    if(test_condition)               % for specular condition 
        post(nk,2) = -post(nk,2);
        post(nk,4) = -post(nk,4);
        
    else 
        
        
        post(nk,2) = 0;                               %Diffusive condition has been met
        Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
        
        theta = rand([sum(nk),1])*2*pi;
        post(nk,3) = Z.*cos(theta);
        post(nk,4) = abs(Z.*sin(theta));
        
        
        
    end

    
   
    
    
  %Figuring out if the particles have moved into to the box. Updates the
  %position and restores the location of the particles
    
    for nk=1:electron_Num
        bottle_neck = box_no(post(nk,1:2), boxes);
        
        
        % Checking for the collision with a box and determining the
        % location of the box of collision. 
        
        while(bottle_neck ~= 0)
            
            dist_X = 0;                  %Finding and updating the X position
            
            X_updated = 0;
            
            
            if(post(nk,3) > 0)
                
                dist_X = post(nk,1) - boxes(bottle_neck,1);
                X_updated = boxes(bottle_neck,1);
                
            else
                
                dist_X = boxes(bottle_neck,2) - post(nk,1);
                X_updated = boxes(bottle_neck,2);
                
                
            end

            dist_Y = 0;                  %Finding and updating the Y position
            Y_updated = 0;
            
            if(post(nk,4) > 0)
                
                dist_Y = post(nk,2) - boxes(bottle_neck, 3);
                Y_updated = boxes(bottle_neck, 3);
                
            else
                
                dist_Y = boxes(bottle_neck, 4) - post(nk,2);
                Y_updated = boxes(bottle_neck, 4);
                
            end

            if(dist_X < dist_Y)
                
                post(nk,1) = X_updated;
                
                if(~spec_boxes(bottle_neck))
                    
                    sgn = -sign(post(nk,3));
                    Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
                    
                    theta = rand()*2*pi;
                    post(nk,3) = sgn.*abs(Z.*cos(theta));
                    post(nk,4) = Z.*sin(theta);
                    
                    
                else 
                    
                    %For specular condition
                    
                    post(nk,3) = -post(nk,3);
                    
                end
                
                
            else
                
                
                post(nk,2) = Y_updated;
                if(~spec_boxes(bottle_neck))
                    
                    sgn = -sign(post(nk,4));
                    Z = sqrt(post(nk,3).^2 + post(nk,4).^2);
                    theta = rand()*2*pi;
                    
                    post(nk,3) = Z.*cos(theta);
                    post(nk,4) = sgn.*abs(Z.*sin(theta));
                    
                else 
                    
                    %For speuclar condition
                    
                    post(nk,4) = -post(nk,4);
                    
                end
            end
            

            bottle_neck = box_no(post(nk,1:2), boxes);
            
            
        end
        
    end


    
    nk = rand(electron_Population, 1) < P_scat;    %Scattering particles using the exponential scattering probability 
    post(nk,3:4) = random(Distr_MB, [sum(nk),2]);

    
    %Calculating and Recording the Temperature of the electrons
   Temp(ni) = (sum(post(:,3).^2) + sum(post(:,4).^2))*Cm/K/2/electron_Population;

    
   %Plotting the electron trajectories and temeperature variance
    for nk=1:electron_Num
       Trj(ni, (2*nk):(2*nk+1)) = post(nk, 1:2);
    end 

    % The animation updates for every 10 iterations
    if animation_plot && mod(ni,10) == 0
        
        figure(6);
        plot(3,1,1);
        hold off;
        plot(post(1:electron_Num,1)./st_size, post(1:electron_Num,2)./st_size, 'o');
        hold on;

        
        
        for nk=1:size(boxes,1)          %Plotting the rectangular boxes 
            
           plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)]./st_size,...
               [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)]./st_size, 'k-');
           
        end

        axis([0 x_R/st_size 0 y_R/st_size]);
        title(sprintf('Electron Trajectories, SM_101092039', electron_Num, electron_Population));
        xlabel('x (nm)');
        ylabel('y (nm)');
        
        
        %Plotting the temeperature variance of the semiconductor 
        if ni > 1
            figure(7)
            plot(3,1,1);
            hold off;
            plot(Tstep*(0:ni-1), Temp(1:ni));
            
             set(gca,'Color', [0 0 0]);
              a_x = gca; 
           
              a_x.GridAlpha = 0.5;  
              a_x.GridColor = [1, 1, 1];
            
            axis([0 Tstep*iter min(Temp(1:ni))*0.98 max(Temp)*1.01]);
            title('Silicon Temperature, SM101092039');
            xlabel('T (s)');
            ylabel('Temperature (K)');
            grid on; 
            
            
        end

       pause(0.05);
    end
end

%Final plotting of the trajectories after all the iterations have been
%completed

figure(6);
subplot(3,1,1);
title(sprintf('Electron Trajectories, SM101092039',electron_Num, electron_Population));
xlabel('Xp(nm)');
ylabel('Yp (nm)');
axis([0 x_R/st_size 0 y_R/st_size]);     
hold on;


%Storing the trajectory after completion 

for ni=1:electron_Num
    
   
    plot(Trj(:,ni*2)./st_size, Trj(:,ni*2+1)./st_size, '.');

end




%Final Plotting of the boxes after completion of iterations
for nk=1:size(boxes,1)
   plot([boxes(nk, 1) boxes(nk, 1) boxes(nk, 2) boxes(nk, 2) boxes(nk, 1)]./st_size,...
       [boxes(nk, 3) boxes(nk, 4) boxes(nk, 4) boxes(nk, 3) boxes(nk, 3)]./st_size, 'k-');
end


%Plotting the temperature of the Semiconductor after the final iteration 
if(~animation_plot)
    figure(7)
    subplot(3,1,1);
    hold off;
    plot(Tstep*(0:iter-1), Temp);
    
     set(gca,'Color', [0 0 0]);
     a_x = gca; 
           
     a_x.GridAlpha = 0.5;  
     a_x.GridColor = [1, 1, 1];
            
    axis([0 Tstep*iter min(Temp)*0.92 max(Temp)*1.01]);
    title('Silicon Temperature, SM101092039');
    xlabel('T(s)');
    ylabel('Temperature (K)');
    grid on; 
    
    
end



density = hist3(post(:,1:2),[200 100])'; %Utilzing the bins from the historgram to plot the the Temperature and Density map 


hist_bins = 10;     %Number of bins utilized
sigma = 3;
[x_elec y_elec]=meshgrid(round(-hist_bins/2):round(hist_bins/2), round(-hist_bins/2):round(hist_bins/2));

M=exp(-x_elec.^2/(2*sigma^2)-y_elec.^2/(2*sigma^2));
M=M./sum(M(:));
figure(8);

imagesc(conv2(density,M,'same'));   %2D convolution for generating a 2D map;
title('Electron Density, SM101092039');

xlabel('Xp (nm)');
ylabel('Yp (nm)');



X_Total_Temp = zeros(ceil(x_R/st_size),ceil(y_R/st_size));
Y_Total_Temp = zeros(ceil(x_R/st_size),ceil(y_R/st_size));                    % Temperature at different locations is being claculated. Velocities are being stored
Total_Temp = zeros(ceil(x_R/st_size),ceil(y_R/st_size));                       %in the bhistorgram bins for calculation



for ni=1:electron_Num             %Figuring out the velocities of the electron 
    
    
  
    
    % Dtermination of the bin 
    x_elec = floor(post(ni,1)/st_size);
    y_elec = floor(post(ni,2)/st_size);
    
    
    if (x_elec==0)
        x_elec = 1;
        
        
    end
    
    if (y_elec==0)
        y_elec= 1;
        
        
    end

    
    Y_Total_Temp(x_elec,y_elec) =  Y_Total_Temp(x_elec,y_elec) + post(ni,3)^2;
    
    X_Total_Temp(x_elec,y_elec) =   X_Total_Temp(x_elec,y_elec) + post(ni,4)^2;  
    
    %Summing all the components of the velocity 
    Total_Temp(x_elec,y_elec)   =      Total_Temp(x_elec,y_elec) + 1;
    
    
    
end




%The temperature can be calculated as:

Temp_sum = (X_Total_Temp + Y_Total_Temp).*Cm./K./2./Total_Temp;
  
Temp_sum = Temp_sum;     


%Plotting the density map

hist_bins = 10;     
sigma = 3;
[x_elec y_elec]=meshgrid(round(-hist_bins/2):round(hist_bins/2), round(-hist_bins/2):round(hist_bins/2)); 

M=exp(-x_elec.^2/(2*sigma^2)-y_elec.^2/(2*(sigma*sigma)));

M=M./sum(M(:));
figure(9);
imagesc(conv2(Temp_sum,M,'same'));       %Performing 2D convolution 


title('Temperature Map, SM101092039');

xlabel('Xp (nm)');
ylabel('Yp (nm)');




function box_num = box_no(pos, boxes)   


%%function for determining the position of the position of the boxes and plotting

    box_num = 0;  %initializing at 0 
    
    for i=1:size(boxes,1)
        
        if(pos(1) > boxes(i,1) && pos(1) < boxes(i,2) && pos(2) > boxes(i,3) && pos(2) < boxes(i,4))
            
            box_num = i;
            
            return;
            
        end
        
    end
    
end
