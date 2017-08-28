% Matlab script ERateDelta18O.m
% Todd Ehlers/Mirjam Schaller  August 19, 2004
%
% This script generates an erosion history from d18O values taking 
% into account different intitial erosion rates and amplitudes and
% writes results to a file for simulation of cosmogenic erosion histories


clear all
close all

% ***
% *** START USER Defined Parameters for Input Erosion Rate
% ***

    % ***
    % *** Define Parameters for Input Erosion Rate Calculation, e.g. d18O value
    % ***

        erate_range = [0.2];            % Define the erosion rates that will be looped through [mm/yr]
        amp_range = [0.5];             % Amplitudes used in calculation to reduce changes in erosion rates
        timestep0 =1000;                   % Time step for resampling and depth calculation
       
        load zachosd18O3Ma.prn;                % Load file with time and d18) values
        t = (zachosd18O3Ma(:,1)')*1000000;     % Assign time values to t in years
        d18O = zachosd18O3Ma(:,2)';            % Assign d18O values to d18O
        
 
    
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
        rho = 2.7;                                    % Eock density

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
% *** Two loops calculating different combinations of amplitude,
% and erosion rates. Save data to files
% ***

count = 1;                               % count amount of columns produced in file

for v = 1: length(amp_range)             % Loop through different amplitudes of interest
    amp = amp_range(v);
    
    amp_out(1) = -666;                    % Assign information about used amplitude to column for file amp_out.out                                    %amplitutde
    amp_out(count+1) = amp;

    
for h = 1:length(erate_range)            % Loop through different erosion rates of interest
    erateinitial = erate_range(h);
    
    erateinitial_out(1) = -666;
    erateinitial_out(count+1) = erateinitial;      % Assign information about used initial erate to column for file erateinitial_out.out
    
    count = count + 2;
    

  % ***
  % *** Calculate Input Erosion Rate and times teps
  % ***
        
        maxd18O = max(d18O);                                 % Search for maximum d18O value
        mind18O = min(d18O);                                 % Search for maximum d18O value
        meand18O = mean(d18O);                               % Determine the mean value of uncorrected d18O values   
        
        normd18O = (d18O - mind18O)/(maxd18O - mind18O);     % Correct for negative values and normalize d18O values to values between 0 and 1   
        meannormd18O = mean(normd18O);                       % Determine the mean value of the normlaized d18O values
        
        % *** 
        % *** Interpolate d18O record to a regular and finely spaced record
        % ***
            % Define new time record
            intert = [0:timestep0:(max(t))];
        
            % Interpolate new d18O record with new time history
            internormd18O = interp1(t,normd18O,intert,'pchip');
        
            
        % *** 
        % *** Transform the interpolated d18O record into a erosion rate
        % *** based on the initial erosion rate and the amplitude
        % ***

        mean_record = mean(amp*internormd18O);
        
        for s = 1:length(intert)
            eratetrue(s) = (erateinitial - mean_record)+(amp*internormd18O(s));
        end
        
       

    % ***
    % *** Calculate Cosmogenic Erosion Rate
    % ***

        m = intert;                           

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
   
        calc_erate (:,1) = intert';                    % Assign Time t to column 1
        calc_erate (:,count-1) = eratetrue';       % Assign true erosion rate eratetrue to column 2,4,6,8 etc.
        calc_erate (:,count) = erosrate';          % Assign cosmogenic erosion rate erosrate to column 3,5,7,9 etc.
   
        
  % ***
  % *** Plot Results of Cosmogenic Erosion Rate Functions
  % *** 
    figure
    
        subplot(211)
        plot(intert,eratetrue,'r-')
        hold on
        plot(intert,erosrate,'b-.')
        title(['Erosion Histories based on d180 values: Amplitude = ' num2str(amp)], 'HorizontalAlignment','center')
        text(((0.8*max(intert))),(1.8*amp),['E Rate = ' num2str(erate_range(1)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        %xlabel('Time [yr]')
        ylabel('Erosion Rate [mm/yr]')
        axis([0 max(intert) 0 (amp_range(1))])
        grid on
        %legend('True Erosion Rate','Cosmogenic Erosion Rate')
    
        subplot(212)
        plot(intert,(erosrate-eratetrue),'r-')
        hold on
        plot(intert,((erosrate-eratetrue)./eratetrue),'b-.')
        %text((0.8*(max(intert))),(1.8*erateinitial),['E Rate = ' num2str(erate_range(2)) ' mm/yr'],'BackgroundColor',[1.0 1.0 1.0], 'HorizontalAlignment','left')
        xlabel('Time [yr]')
        ylabel('Difference True-Cosmogenic Erosion Rate')
        axis([0 max(intert) (-amp_range(1)/4) (amp_range(1)/4)])
        grid on
        legend('Difference abosolut [mm/yr]','Difference relative')
    

end             %end h loop

end             %end v loop


% ***
% ***END OF PROGRAM TO CALCULATE TRUE AND COSMOGENIC EROSION RATE
% ***

% ***
% ***Save data to different files
% ***

    save calc_erate.out calc_erate -ascii       %save all values of time, eratetrue and erosrate to one file
    
    save amp_out.out amp_out -ascii             % file storing information about used amplitude 
    
    save erateinitial_out.out erateinitial_out -ascii      % file storing information about used initial erosion rate 
   
 
    
           % ***
       % *** Plot input d18O and true erosion rate
       % *** 
            
         figure
         subplot(211)
            plot(t,d18O,'.')
            title('Input d18O Record')
            xlabel('Time (yr)')
            ylabel('d18O')
         subplot(212)
            plot(t,normd18O,'.')
            title('Input Normalized Record - Not Interpolated')
            xlabel('Time (yr)')
            ylabel('Normalized d18O')
            
            
          figure
          subplot(211)
            plot(intert,internormd18O,'.')
            title('Interpolated Record ')
            xlabel('Time (yr)')
            ylabel('Interpolated Normalized d18O')
          subplot(212)
            plot(intert,eratetrue)
            xlabel('Time Before Present (yrs)')
            ylabel('Erosion Rate (mm/yr)')
            title('Erosion History Based on d18O Record [Zachos, 2001]')
    
 