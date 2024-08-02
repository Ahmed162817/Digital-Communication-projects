clc
clear
close all;

%% ------------------------Matched filter & correlator part-------------------------%%

Ts = 1;                        % it refer to the symbol duration Ts = 1 sec
% generation of pulse p[n] but divided over sqrt(55) to normalize it's energy 
pulse = [5 4 3 2 1]/sqrt(55);  
Number_of_bits = 10;           % this variable represent the number of the random bits that are generated
Data = randi([0 1] , 1 , Number_of_bits);    % array contain the random data bits (0 & 1)

Amplitude = 1;         % it refer to the amplitude of the impulse train 
PAM_signal = (2 * Data - 1) * Amplitude;      % This is the PAM sampling which represent (+Amplitude or -Amplitude) every Ts = 1 sec
impulse_train = upsample(PAM_signal,5);       % we upsample the input signal (PAM_signal) by inserting 4 zeros (5-1) in each symbol duration Ts 

Tx_out = conv(impulse_train,pulse);     % where Tx_out refer to the signal at the output of the transmitter 
% where the length(Tx_out) = length(impulse_train) + length(pulse) - 1     [According to the discrete convolution]

% For plotting, we will use the stem() function instead of the plot() function as stem() is specified for the discrete plotting
figure;
subplot(2,1,1);
stem((0:Ts/5:(length(impulse_train)-1)/5), impulse_train, 'r', 'LineWidth', 2);
title('Impulse Train');
xlabel('sampling Time instants');
grid on;

subplot(2,1,2);
plot(0:length(Tx_out)-1, Tx_out, 'b', 'LineWidth', 2);
title('Output signal of the Transmitter (Tx_out)');
xlabel('Time');
grid on;

%------------- Methods of filtering the Tx_out signal (2 Methods)------%

% First Method is by take Tx_out & apply convolution with impulse response of the Matched filter h[n] = Pulse[Ts - n]     (Matched Filter)

impulse_response_MF1 = fliplr(pulse);             % Matched filter impulse response
Rx_out1 = conv(Tx_out, impulse_response_MF1);      % Rx_out is the output signal of the matched filter

% second Method is by take Tx_out & apply convolution with impulse response of another filter h(t) = 1 from 0 --> Ts    (NOT Matched Filter)

impulse_response_MF2 = [1 1 1 1 1]/sqrt(5);       % First step we must normalize the impulse response to obtain unity energy
Rx_out2 = conv(Tx_out, impulse_response_MF2);

% we make subplotting to compare between the two outputs of two different filters (Rx_out1,Rx_out2)

figure;
subplot(2,1,1);
plot(0:length(Rx_out1)-1 , Rx_out1 , 'k', 'LineWidth', 2);        
title('Output of the first Filter that its impulse response h[n] = pulse[Ts-n]');
xlabel('Time');
grid on;

subplot(2,1,2);
plot(0:length(Rx_out2)-1, Rx_out2, 'r', 'LineWidth', 2);
title('Output of the second Filter that its impulse response h[n] = 1 from 0 ---> Ts');
xlabel('Time');
grid on;

figure;
subplot(2,1,1);
stem(0:length(Rx_out1)-1 , Rx_out1 , 'k', 'LineWidth', 2);
title('Output of The First Matched Filter in Discrete Form');
xlabel('Time');
grid on;

subplot(2,1,2);
stem(0:length(Rx_out2)-1, Rx_out2, 'r', 'LineWidth', 2);
title('Output of the second Filter in Discrete Form');
xlabel('Time');
grid on;

% Build the correlator & compare it's output with the output of the Matched Filter

correlator_out = zeros(1,length(Number_of_bits));  %the output signal of the correlator is the same length as the number of random bits generated 

j = 1;            % j variable is used to increment the index of correlator_out
for i = 1 : 5 : length(Tx_out)-5
        % Calculate the correlation for each segment (5 samples) of Tx_out and update correlator_out
        correlator_out(1,j) = sum( pulse .* Tx_out(1,i:i+4) );
        j = j + 1;
end

% plot the output of the correlator and compare it with the Matched Filter
figure;
subplot(2,2,[1 2]);
plot(0:length(Rx_out1)-1 , Rx_out1 , 'g', 'LineWidth', 2);
title('Output of The Matched Filter');
xlabel('Time');
grid on;

subplot(2,2,3);
plot(0:length(correlator_out)-1, correlator_out , 'r', 'LineWidth', 2);
title('Output of the Correlator');
xlabel('Time');
grid on;

subplot(2,2,4);
stem(0:length(correlator_out)-1, correlator_out , 'b', 'LineWidth', 2);
title('Output of the Correlator in Discrete Form');
xlabel('Time');
grid on;

%% ---------------------------------Noise Analysis part-------------------------- %%
%  a):
Number_of_bits = 10000;           % this variable represent the number of the random bits that are generated
pulse = [5 4 3 2 1]/sqrt(55);  
Data = randi([0 1] , 1 , Number_of_bits);    % array contain the random data bits (0 & 1)
Amplitude = 1;         % it refer to the amplitude of the impulse train 
PAM_signal = (2 * Data - 1) * Amplitude;      % This is the PAM sampling which represent (+Amplitude or -Amplitude) every Ts = 1 sec
impulse_train = upsample(PAM_signal,5);       % we upsample the input signal (PAM_signal) by inserting 4 zeros (5-1) in each symbol duration Ts 
Tx_out = conv(impulse_train,pulse);     % where Tx_out refer to the signal at the output of the transmitter 
impulse_response_MF1 = fliplr(pulse);             % Matched filter impulse response
impulse_response_MF2 = [5 5 5 5 5]/sqrt(125);       % First step we must normalize the impulse response to obtain unity energy

% b):
Gaussian_Noise = randn( 1,length(Tx_out) );  % generate gaussian white noise of zero mean and unity variance of length equal to that of transmitted signal

% c):
Noise_Variance = [1.585 1.25 1 0.794 0.631 0.501 0.3981 0.316 0.2511 0.1995 0.15845]  ;  % Variable to control gaussian noise variance value

% d):
BER_Simulated = zeros(2,length(Noise_Variance));
BER_Theoretical = zeros(1,length(Noise_Variance));
for j= 1 :length(Noise_Variance)
    Noisy_Tx_out = Tx_out + ( Gaussian_Noise * sqrt( Noise_Variance(j)/2 )) ;
    Rx_out1 = conv(Tx_out, impulse_response_MF1);      % Noisy_Rx_out1 is the output signal of the matched filter with noise added
    Noisy_Rx_out1 = conv(Noisy_Tx_out, impulse_response_MF1);      % Noisy_Rx_out1 is the output signal of the matched filter with noise added
    Noisy_Rx_out2 = conv(Noisy_Tx_out, impulse_response_MF2);
    Sampled_Output1 = zeros(1,Number_of_bits); %Samples Noisy_Rx_out1
    Sampled_Output2 = zeros(1,Number_of_bits); %Samples Noisy_Rx_out2

    Sampled_Output_Index = 1;
    for i= 5:5: (Number_of_bits*5)
        Sampled_Output1(Sampled_Output_Index ) = Noisy_Rx_out1(i);
        Sampled_Output2(Sampled_Output_Index ) = Noisy_Rx_out2(i);
        Sampled_Output_Index = Sampled_Output_Index + 1;
    end
    Number_Of_Errors = zeros(1,2);
    %Calculation Of BER
    for i = 1:Number_of_bits   
        if(PAM_signal(i) == 1 && Sampled_Output1(i) <=0 )
            Number_Of_Errors(1) = Number_Of_Errors(1) + 1;
        elseif(PAM_signal(i) == -1 && Sampled_Output1(i) >0 )
            Number_Of_Errors(1) = Number_Of_Errors(1) + 1;
        end 
        if(PAM_signal(i) == 1 && Sampled_Output2(i) <=0 )
            Number_Of_Errors(2) = Number_Of_Errors(2) + 1;
        elseif(PAM_signal(i) == -1 && Sampled_Output2(i) >0 )
            Number_Of_Errors(2) = Number_Of_Errors(2) + 1;
        end 
    end
    BER_Simulated(1,j) = Number_Of_Errors(1) / Number_of_bits  ;
    BER_Simulated(2,j) = Number_Of_Errors(2) / Number_of_bits  ;
    BER_Theoretical(1,j) = 0.5 * erfc(sqrt(1/Noise_Variance(j) ) ); 
end
  
%plotting BER vs SNR
figure;
SNR_DB=-2:length(Noise_Variance)-3;
semilogy(SNR_DB, BER_Simulated(1,:) ,'r', 'linewidth', 2)
hold on;
semilogy(SNR_DB, BER_Simulated(2,:),'b', 'linewidth', 2 )
hold on;
semilogy(SNR_DB, BER_Theoretical , 'g', 'linewidth', 2)
legend('Matched Filter','Ideal fitler','Theoretical Value')
grid on
xlabel('SNR Value in dB')
ylabel('BER Value')
title('Theoretical and Matched Filter and Ideal Filter BER vs SNR in dB ')

%% ***************** ISI and Raised cosine filter part*****************%%
% Parameters
numBits = 100;          % Number of bits
rolloff = [0,0,1,1];          % Roll-off factor
delay=[2,8,2,8];            %Delay introduced to filter coefficient
sps =5;       %Samples taken per symbol preferred large to increase resolution
Fs = 1;     %Sampling freq  where each bit is represented by a symbol

% Generate random bits
dataTx = randi([0,1],1,numBits)*2-1;  % Randomize bits to -1 or 1
dataTx_upsampled = upsample(dataTx,sps);

for i=1:4
    
    % Transmitter
        %Generating the filter coefficients 
            srrc = rcosine(Fs ,sps*Fs, 'sqrt' , rolloff(i), delay(i) );
        
       % Apply SRRC filter to the upsampled data by convolution
            txSig = conv(dataTx_upsampled, srrc);
        
        %Calculating the start of the main intended output     
        start = length(srrc)+1;
        
    % Eyediagram at A(transmitter)  
        eyediagram(txSig(start-1:end-start-sps-1),2*sps);
        xlabel('Time (samples)');
        title(['Filter output at A with rolloff and delay equal  ' num2str(rolloff(i)) ' & ' num2str(delay(i))]);
        
    % Receiver
        %Using the same filter srrc as the transmitter
            RxSig = conv(txSig,srrc);  
                 
        % Eyediagram at B(reciever)
           eyediagram(RxSig(start-1:end-start-sps-1),2*sps);
           xlabel('Time (samples)');
           title(['Filter output at B with rolloff and delay equal  ' num2str(rolloff(i)) ' & ' num2str(delay(i))]);

end
