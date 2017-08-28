% Matlab script erate_periodicV2.m
% Todd Ehlers/Mirjam Schaller  August 19, 2004
%
% This script generates a sinusoidal erosion history with different
% wavelength and amplitude and writes results to a file for simulation
% of cosmogenic exposure ages and basin thermal histories.
%

clear all
close all

% ***
% *** START USER Defined Parameters for Input Erosion Rate
% ***

    % ***
    % *** Define Parameters for Input Erosion Rate Calculation
    % ***
    
        erate_range=[0.02];       % Define the erosion rates that will be looped through [mm/yr]
        amp_range = [0.1];              % Amplitudes of periodic function
        wave_range = [41000 100000];        % Wavelength of periodic function [yr]
        base = 1.00000000001; 	            % Background, or baseline erosion rate that oscillation occurs around [mm/yr]

        timestep0 =100;                    % time interval at which calculations will be done
        timeend = 7e6;                      % tune in yrs to which calculation should run
        t = 0:timestep0:timeend;            % time span overwhich erosion rates are calculated [yrs]

    
    % ***
    % *** Define Parameters for Cosmogenic Erosion Rate Calculation
    % ***

        prodratenucsurf = 4.431;                 % Nucleonic Production Rate at SLHL
        prodratemuonstoppedsurf = 0.025;        % Stipped Muonic Production Rate at SLHL
        prodratemuonfastsurf = 0.037;           % Fast Muonic Production Rate at SLHL


        a(1) = 0 * prodratenucsurf;        % Production parameter
        a(2) = 1 * prodratenucsurf;
        a(3) = 0.845 * prodratemuonstoppedsurf;
        a(4) = -0.05 * prodratemuonstoppedsurf;
        a(5) = 0.205 * prodratemuonstoppedsurf;
        a(6) = 0.01 * prodratemuonfastsurf;
        a(7) = 0.615 * prodratemuonfastsurf;
        a(8) = 0.375 * prodratemuonfastsurf;

        absorp(1) = 1;                    % Absorption length in g/cm2
        absorp(2) = 157;
        absorp(3) = 1030;
        absorp(4) = 160;
        absorp(5) = 3000;
        absorp(6) = 100;
        absorp(7) = 1520;
        absorp(8) = 7600;

        decay = 4.99746e-07;                           % Decay Constant
        rho = 2.7;                                  % Eock density

        timestep = 0;
        cutoff = 0.001;                         % 1 permille Cutoff for Production Calculation at Depth
        interfacedepth = 46;
    
% ***
% *** END OF USER DEFINED PARAMETERS
% ***
 



% ***
% *** START OF PROGRAM TO CALCULATE TRUE AND COSMOGENIC EROSION RATE
% ***

% ***
% *** Three loops calculating different combinations of ampltitude,
% wave length and erosion rates. Save data to files
% ***

count = 1;                               % count amount of columns produced in file

for v = 1: length(amp_range)             % Loop through different amplitudes of interest
    amp = amp_range(v);
    
    amp_out(1) = -666;                    % Assign information about used amplitude to column for file amp_out.out                                    %amplitutde
    amp_out(count+1) = amp;
    
for w = 1: length(wave_range)            % Loop through different wave length of interest
    wave = wave_range(w);
    
    figure                               % Make new figure
    
    wave_out(1) = -666;                  % Assign information about used wave length to column for file wave_out.out                                    %amplitutde
    wave_out(count+1) = wave;

    
for h = 1:length(erate_range)            % Loop through different erosion rates of interest
    erateinitial = erate_range(h);
    
    erateinitial_out(1) = -666;
    erateinitial_out(count+1) = erateinitial;      % Assign information about used initial erate to column for file erateinitial_out.out
    
    count = count + 2;
    
    
    % ***
    % *** Calculate Input Erosion Rate Function
    % ***

        for i=1:length(t)  
            if i<length(t)/2
                eratetrue(i) =  erateinitial * (base + amp * sin( ((2*pi)/wave * t(i)) ));
            else
             eratetrue(i) =  erateinitial;   
            end
        end

    

    % ***
    % *** Calculate Cosmogenic Erosion Rate
    % ***

        m = t;                           

    % Integrate over True Erosion Rate * dt to calculate transient Depth rock was buried_
    % at time t and calculate Production Rate(Depth), and integrated Nuclide Inventory

             
    for j = 1:length(m)                     
                                        
                prodrateinstantaneous = 0;
                depth = 0;                          % in cm
                concsurf = 0;
            
        for k = j:length(m)
   
                % Subtract one half Step-Depth from current Step(k) so that
                % Calculation is done in the middle of depth intervall and
                % calculate ProdRate at Depth cm
            depth = depth + (eratetrue(k))/10 * timestep0 * 0.5; 

                % Calculate Production Rate at given Depth   
            prodratenuc = a(1) * exp(-rho / absorp(1) * depth) + a(2) * exp(-rho / absorp(2) * depth);
            prodratemuonstopped = a(3) * exp(-rho / absorp(3) * depth) + a(4) * exp(-rho / absorp(4) * depth) + a(5) * exp(-rho / absorp(5) * depth);
            prodratemuonfast = a(6) * exp(-rho / absorp(6) * depth) + a(7) * exp(-rho / absorp(7) * depth) + a(8) * exp(-rho / absorp(8) * depth);

                % Calculate Total Production Rate
                prodrateinstantaneous =prodratenuc+prodratemuonstopped+prodratemuonfast;
       
                % Terminate Loop if Production Rate at Depth < than "Cutoff" *  Prod(0)
                % This saves considerable calculation time for long runs
           % if ((cutoff > 0) & (prodrateinstantaneous < cutoff * (prodratenucsurf + prodratemuonfastsurf + prodratemuonstoppedsurf)));
           if ((cutoff > 0) && (prodrateinstantaneous < cutoff * (prodratenucsurf + prodratemuonfastsurf + prodratemuonstoppedsurf)));
                %display('*******************************************************')
                %display('*******************************************************')
                k;
                break
            end

                % Then Calculate the integrated Nuclide Concentration rock has obtained en route to
                % the Surface where it is eroded concnentration produced in depth corrected for decay 
            concdepth = prodrateinstantaneous * timestep0 * exp(-timestep0 * decay);

                 % Correction for decay during transport
            concdepth = concdepth * exp(-decay * depth / (eratetrue(k)/10));
            
                 % Concentration at surface
            concsurf = concsurf + concdepth;
            
                % Re-Add lost half Step for total depth of subsequent steps
            depth = depth + (eratetrue(k))/10 * timestep0 * 0.5;
           
        
        end                     % end k loop
        
                                    
                % Iterative Calculation of Apparent Erosion Rates
            erosrate(j) = 0.1;       % First Guess for Erosion Rate in cm/yr
            conc = 0;                % First Guess for Concentration
            for q = 1:8             
                conc = conc + (a(q) / (decay + (erosrate(j)/10 * rho / (absorp(q)))));
            end
            
        while abs(conc / concsurf - 1) > 0.0001
            erosrate(j) = erosrate(j) - erosrate(j) * (concsurf - conc) / concsurf;
            
            %Loop for 8 Production Mode
            conc = 0;
            for p = 1:8             
                conc = conc + (a(p) / (decay + (erosrate(j)/10 * rho / (absorp(p)))));
            end
            
        end

  end                       % end j loop
  
  % ***
  % *** Attribute data to different columns
  % *** 
   
        calc_erate (:,1) = t';                    % Assign Time t to column 1
        calc_erate (:,count-1) = eratetrue';       % Assign true erosion rate eratetrue to column 2,4,6,8 etc.
        calc_erate (:,count) = erosrate';          % Assign cosmogenic erosion rate erosrate to column 3,5,7,9 etc.
        
    
  % ***
  % *** Plot Results of Cosmogenic Erosion Rate Functions
  % *** 
    
    if h==1
        subplot(221)
        plot(t,eratetrue,'r-')
        hold on
        plot(t,erosrate,'b-.')
        title(['Erosion Histories: Amplitude = ' num2str(amp) '; Wave Length = ' num2str(wave) ' years'], 'HorizontalAlignment','left')
        text((2.5*wave),(1.8*erateinitial),['E Rate = ' num2str(erate_range(1)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        %xlabel('Time [yr]')
        ylabel('Erosion Rate [mm/yr]')
        axis([0 (4*wave) 0 (2*erateinitial)])
        grid on
        %legend('True Erosion Rate','Cosmogenic Erosion Rate')
    elseif h==2
        subplot(222)
        plot(t,eratetrue,'r-')
        hold on
        plot(t,erosrate,'b-.')
        text((2.5*wave),(1.8*erateinitial),['E Rate = ' num2str(erate_range(2)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        %xlabel('Time [yr]')
        %ylabel('Erosion Rate [mm/yr]')
        axis([0 (4*wave) 0 (2*erateinitial)])
        grid on
        %legend('True Erosion Rate','Cosmogenic Erosion Rate')
    elseif h==3
        subplot(223)
        plot(t,eratetrue,'r-')
        hold on
        plot(t,erosrate,'b-.')
        text((2.5*wave),(1.8*erateinitial),['E Rate = ' num2str(erate_range(3)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        xlabel('Time [yr]')
        ylabel('Erosion Rate [mm/yr]')
        axis([0 (4*wave) 0 (2*erateinitial)])
        grid on
        %legend('True Erosion Rate','Cosmogenic Erosion Rate')
    elseif h==4
        subplot(224)
        plot(t,eratetrue,'r-')
        hold on
        plot(t,erosrate,'b-.')
        text((2.5*wave),(1.8*erateinitial),['E Rate = ' num2str(erate_range(4)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        xlabel('Time [yr]')
        %ylabel('Erosion Rate ]mm/yr]')
        axis([0 (4*wave) 0 (2*erateinitial)])
        grid on
        %legend('True Erosion Rate','Cosmogenic Erosion Rate')
    end         % end if h statement for subplotting output

end             %end h loop

end             %end w loop

end             %end v loop


% ***
% ***END OF PROGRAM TO CALCULATE TRUE AND COSMOGENIC EROSION RATE
% ***

% ***
% ***Save data to different files
% ***

    save calc_erate.out calc_erate -ascii       %save all values of time, eratetrue and erosrate to one file
    
    save amp_out.out amp_out -ascii             % file storing information about used amplitude 
  
    save wave_out.out wave_out -ascii           % file storing information about used wave length 
    
    save erateinitial_out.out erateinitial_out -ascii      % file storing information about used initial erosion rate 

    %save all_variables.mat          % binary matlab file that can be loaded 
                                    % using the load command for plotting