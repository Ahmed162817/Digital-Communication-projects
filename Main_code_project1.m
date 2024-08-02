clc
clear 
close all ;

%% Generating the ensemble with different line codes
A = 4 ;        % Amplitude of the pulse
Number_of_realizations = 500;
Number_of_bits = 100;             % in each realization there will be 100 bits
Number_of_samples = 7;            % there are 7 samples for each bit duration
% since the time of DAC =10ms and Tb (period of each bit) =70ms so, we need 7 samples in each bit duration (Tb)

Data = randi([0 1] , Number_of_realizations , Number_of_bits+1);   % generate all the random bits of the ensemble 
% where we add additional bit to start from ( Time_Delay ---> Number_of_bits+Time_Delay )

Tx1 = Data * A ;           % unipolar line coding
Tx2 = (2*Data-1) * A;     % polar line coding (General for NRZ & RZ)

Tx1_unipolar = repelem(Tx1 , 1 , Number_of_samples);    % repelem() function repeat each element in the array Tx1 7 times [unipolar output signalling]
Tx2_polar_NRZ = repelem(Tx2 , 1 , Number_of_samples);    % Tx2_out represent polar NRZ only 
Tx3_polar_RZ = Tx2_polar_NRZ;                      % initially equate polar return to zero with polar non-return to zero

for i = 1 : Number_of_realizations
  for j = 1 : Number_of_bits+1
     Tx3_polar_RZ( i , j*Number_of_samples - 2 : j*Number_of_samples) = 0;     % assume number of zeros in each bit duration = 3  
  end  
end
   
Time_Delay = randi([0 Number_of_samples-1] , Number_of_realizations , 1);         % Time_Delay refer to the time delay at the start of each realization   

% initialize each line coding with zeros
Tout_unipolar = zeros(Number_of_realizations,Number_of_bits*Number_of_samples);
Tout_polar_NRZ = zeros(Number_of_realizations,Number_of_bits*Number_of_samples);
Tout_polar_RZ = zeros(Number_of_realizations,Number_of_bits*Number_of_samples);
% Adding the time delay (time shift Td) to each line coding
for i = 1 : Number_of_realizations
   Tout_unipolar(i,:) = Tx1_unipolar(i , (Time_Delay(i)+1) : (Number_of_bits*Number_of_samples+Time_Delay(i)));
   Tout_polar_NRZ(i,:) = Tx2_polar_NRZ(i , (Time_Delay(i)+1) : (Number_of_bits*Number_of_samples+Time_Delay(i)));
   Tout_polar_RZ(i,:) = Tx3_polar_RZ(i , (Time_Delay(i)+1) : (Number_of_bits*Number_of_samples+Time_Delay(i)));
end
  
%% calculating the Statistical Mean for different line codes

% for unipolar line code
Statistical_mean_unipolar = zeros(1,Number_of_bits*Number_of_samples);
for j = 1 : (Number_of_bits*Number_of_samples)
   Statistical_mean_unipolar(1,j) = Statistical_mean_unipolar(1,j) +  (sum(Tout_unipolar(:,j))/Number_of_realizations); 
end
 % Plotting the statistical mean of unipolar line code
 time = 1 : (Number_of_bits * Number_of_samples) ;   % time vector refer to the x-axis 
    figure;
    subplot(3,1,1);
    plot( time , Statistical_mean_unipolar,'r')
    xlabel('Time');
    title('Statisitcal Mean of Unipolar Line code');
   
 % for polar NRZ line code
Statistical_mean_polar_NRZ = zeros(1,Number_of_bits*Number_of_samples);
for j = 1 : (Number_of_bits*Number_of_samples)
   Statistical_mean_polar_NRZ(1,j) = Statistical_mean_polar_NRZ(1,j) +  (sum(Tout_polar_NRZ(:,j))/Number_of_realizations); 
end
 % Plotting the statistical mean of polar NRZ line code 
    subplot(3,1,2);
    plot( time , Statistical_mean_polar_NRZ,'b')
    xlabel('Time');
    title('Statisitcal Mean of Polar NRZ Line code');  
 
 % for polar RZ line code
Statistical_mean_polar_RZ = zeros(1,Number_of_bits*Number_of_samples);
for j = 1 : (Number_of_bits*Number_of_samples)
   Statistical_mean_polar_RZ(1,j) = Statistical_mean_polar_RZ(1,j) + (sum(Tout_polar_RZ(:,j))/Number_of_realizations); 
end
 % Plotting the statistical mean of polar NRZ line code 
    subplot(3,1,3);
    plot( time , Statistical_mean_polar_NRZ,'k')
    xlabel('Time');
    title('Statisitcal Mean of Polar RZ Line code'); 
    
%% calculating the Time Mean for different line codes

%for unipolar line code
time_mean_unipolar = zeros(Number_of_realizations,1);
for i = 1 : Number_of_realizations
  time_mean_unipolar(i,1) =  time_mean_unipolar(i,1) + sum(Tout_unipolar(i,:))/(Number_of_bits*Number_of_samples); 
end

%for polar NRZ line code
time_mean_polar_NRZ = zeros(Number_of_realizations,1);
for i = 1 : Number_of_realizations
  time_mean_polar_NRZ(i,1) =  time_mean_polar_NRZ(i,1) + sum(Tout_polar_NRZ(i,:))/(Number_of_bits*Number_of_samples); 
end

%for polar RZ line code
time_mean_polar_RZ = zeros(Number_of_realizations,1);
for i = 1 : Number_of_realizations
  time_mean_polar_RZ(i,1) =  time_mean_polar_RZ(i,1) + sum(Tout_polar_RZ(i,:))/(Number_of_bits*Number_of_samples); 
end

% Plotting the time mean of unipolar line code
Realization_number = 1 : Number_of_realizations ; 
    figure;
    subplot(3,1,1);
    plot( Realization_number , time_mean_unipolar,'c')
    xlabel('Number of Realizations');
    title('Time Mean of Unipolar Line code');

% Plotting the time mean of polar NRZ line code      
    subplot(3,1,2);
    plot( Realization_number , time_mean_polar_NRZ ,'y')
    xlabel('Number of Realizations');
    title('Time Mean of polar NRZ Line code');


% Plotting the time mean of polar RZ line code      
    subplot(3,1,3);
    plot( Realization_number , time_mean_polar_RZ ,'r')
    xlabel('Number of Realizations');
    title('Time Mean of polar RZ Line code');

%% calculating the Time Autocorrelation Function for different line codes  (For single realization only)

% while(1)
%     realization_value = input('please enter the realization number range from 1 --> 500to bring its time autocorrelation with different line codes : ');
%     if (realization_value < 1 || realization_value > 500)
%         fprintf('invalid value , please try again');
%     else
%         break;
%     end
% end

% for unipolar line code
time_ACF_unipolar = time_ACF(Tout_unipolar,1);        % calling the autocorrelation function for different line codes
time_ACF_polar_NRZ = time_ACF(Tout_polar_NRZ,1);     
time_ACF_polar_RZ = time_ACF(Tout_polar_RZ,1);

% Plotting the time Autocorrelation of unipolar line code
time_range = -(Number_of_bits*Number_of_samples)/2 : (Number_of_bits*Number_of_samples-1)/2;
    figure;
    subplot(3,1,1);
    plot(time_range , time_ACF_unipolar, 'k');   % we concatenate the time autocorrelation after flipping it as the ACF is even function
    xlabel('Time');
    title('Time Autocorrelation of unipolar');

% Plotting the time Autocorrelation of polar NRZ line code
    subplot(3,1,2);
    plot(time_range , time_ACF_polar_NRZ , 'g');  
    xlabel('Time');
    title('Time Autocorrelation of polar NRZ');

% Plotting the time Autocorrelation of polar RZ line code
    subplot(3,1,3);
    plot(time_range , time_ACF_polar_RZ , 'b');  
    xlabel('Time');
    title('Time Autocorrelation of polar RZ');

%% calculating the Statistical Autocorrelation Function for different line codes

stat_ACF_unipolar = stat_ACF(Tout_unipolar);       
stat_ACF_polar_NRZ = stat_ACF(Tout_polar_NRZ);        
stat_ACF_polar_RZ = stat_ACF(Tout_polar_RZ);

% Plotting the statistical Autocorrelation of unipolar line code
    figure;
    subplot(3,1,1);
    plot(time_range , stat_ACF_unipolar , 'k');  
    xlabel('Time');
    title('Statistical Autocorrelation of Unipolar');
 
% Plotting the statistical Autocorrelation of polar NRZ line code
    subplot(3,1,2);
    plot(time_range , stat_ACF_polar_NRZ , 'r');  
    xlabel('Time');
    title('Statistical Autocorrelation of polar NRZ');

% Plotting the statistical Autocorrelation of polar RZ line code
    subplot(3,1,3);
    plot(time_range , stat_ACF_polar_RZ , 'g');  
    xlabel('Time');
    title('Statistical Autocorrelation of polar RZ');
       
%% calculating the power spectral density (PSD) & Bandwidth for each line code

Ts = 10 * 10^(-3);        % the sampling period obtained from the DAC which is 10 ms
Fs = 1/Ts;                % the sampling freq = 100 Hz 
Length_ACF = length(stat_ACF_unipolar);   % Length_ACF is general length for unipolar & polar NRZ & polar RZ

% For unipolar line code

PSD_unipolar = fft(stat_ACF_unipolar,Length_ACF);

    freq_range = (-Length_ACF/2 : Length_ACF/2-1) * (Fs/Length_ACF);
    figure;
    subplot(3,1,1);
    plot(freq_range , abs(fftshift(PSD_unipolar))/Length_ACF , 'r');  % we use the fft shift to make the PSD symmetric around zero
    xlabel('Frequency in Hz');
    title('PSD of the unipolar line code');    
 
% For polar NRZ line code
PSD_polar_NRZ = fft(stat_ACF_polar_NRZ,Length_ACF);
    subplot(3,1,2);
    plot(freq_range , abs(fftshift(PSD_polar_NRZ))/Length_ACF , 'b');  
    xlabel('Frequency in Hz');
    title('PSD of the polar NRZ line code');

% For polar RZ line code
PSD_polar_RZ = fft(stat_ACF_polar_RZ,Length_ACF);
    subplot(3,1,3);
    plot(freq_range , abs(fftshift(PSD_polar_RZ))/Length_ACF , 'k');  
    xlabel('Frequency in Hz');
    title('PSD of the polar RZ line code'); 
    
%% STATISTICAL AUTOCORRELATION FUNCTION IMPLEMENTATION
function stat_autocorr = stat_ACF(ensemble_data)    
% ensemble_data refer to the input line code data (unipolar, polarNRZ, polarRZ) 
[m, n] = size(ensemble_data);   
% where the m represent the number of realizations & n represent the (no. of bits * no. of samples)
stat_autocorr = zeros(1,n);

    for time_shift = 0 : (n-1)                
       shifted_data = circshift(ensemble_data,time_shift+n/2,2);      
       % This function shift the ensemble_data by time_shift positions and
       %this shift is for 2nd dimension only (applied for columns)
       stat_autocorr(1,time_shift+1) = sum(sum(ensemble_data .* shifted_data)/n) / m; 
    end
    
end

%% TIME AUTOCORRELATION FUNCTION IMPLEMENTATION
function time_autocorr = time_ACF(ensemble_data,realization_number)    
% realization number is an input argument to specify the desired
%realization for time autocorrelation

Number_of_columns = size(ensemble_data,2);
    time_autocorr = zeros(1,Number_of_columns);
%Starting shift from mid position of realization to get
% Positive & negative time differences 
Target_Realization = ensemble_data(realization_number,:); 

  for time_diff = 0 : (Number_of_columns-1)
      
       % rlz_shifted represent the shifted version of the realization
       rlz_shifted = circshift(Target_Realization,[0,time_diff-Number_of_columns/2]);   
      
       product = Target_Realization .* rlz_shifted ;
       time_autocorr(1,time_diff+1) = sum(product);
  end
time_autocorr = time_autocorr / Number_of_columns ;
end
 