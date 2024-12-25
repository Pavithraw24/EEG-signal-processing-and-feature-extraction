
warning('off', 'EDFANNOT2EVT:M');  
warning('off', 'SOPEN:OverflowDetection');  


addpath('E:\dss');
savepath;

% Specify the paths to your EDF files
healthy_file = 'E:\dss\H S9 EC.edf';
depressed_file = 'E:\dss\MDD S9 EC.edf';


[hdr_healthy, record_healthy] = sopen(healthy_file);
[hdr_depressed, record_depressed] = sopen(depressed_file);


signal_healthy = sread(hdr_healthy);  % Extract EEG signal for healthy data
signal_depressed = sread(hdr_depressed);  % Extract EEG signal for depressed data

sclose(hdr_healthy);
sclose(hdr_depressed);

disp('Size of Healthy EEG Signal:');
disp(size(signal_healthy));  

disp('Size of Depressed EEG Signal:');
disp(size(signal_depressed));

num_channels_healthy = size(signal_healthy, 2);  
figure;
tiledlayout(ceil(num_channels_healthy / 2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:num_channels_healthy
    nexttile;
    plot(signal_healthy(:, i));  
    title(['Healthy EEG - Ch ', num2str(i)]);
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Healthy EEG Signals');


num_channels_depressed = size(signal_depressed, 2);  
figure;
tiledlayout(ceil(num_channels_depressed / 2), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:num_channels_depressed
    nexttile;
    plot(signal_depressed(:, i));  
    title(['Depressed EEG - Ch ', num2str(i)]);
    xlabel('Samples');
    ylabel('Amplitude');
    grid on;
end
sgtitle('Depressed EEG Signals'); 

sampling_rate_healthy = hdr_healthy.SampleRate;
sampling_rate_depressed = hdr_depressed.SampleRate;  


disp(['Sampling Rate for Healthy Data: ', num2str(sampling_rate_healthy), ' Hz']);
disp(['Sampling Rate for Depressed Data: ', num2str(sampling_rate_depressed), ' Hz']);

bands = [0.5 4; 4 8; 8 13; 13 30];  
num_bands = size(bands, 1);  

% Preallocate features matrices correctly
features_healthy = zeros(num_channels_healthy, num_bands);
features_depressed = zeros(num_channels_depressed, num_bands);

% Calculate band powers for healthy signals
for i = 1:num_channels_healthy
    disp(['Processing Healthy Channel: ', num2str(i)]);
    
    
    [pxx_healthy, f_healthy] = pwelch(signal_healthy(:, i), [], [], [], 256); 
    for j = 1:num_bands
        band_freqs = bands(j, :);  
        
       
        idx_band = find(f_healthy >= band_freqs(1) & f_healthy <= band_freqs(2));
        
        
        if ~isempty(idx_band)
            features_healthy(i, j) = trapz(f_healthy(idx_band), pxx_healthy(idx_band));
        else
            features_healthy(i, j) = 0;
        end
    end
end


for i = 1:num_channels_depressed
    disp(['Processing Depressed Channel: ', num2str(i)]);
   
    [pxx_depressed, f_depressed] = pwelch(signal_depressed(:, i), [], [], [], 256);   

   
    for j = 1:num_bands
        band_freqs = bands(j, :);  
        
     
        idx_band = find(f_depressed >= band_freqs(1) & f_depressed <= band_freqs(2));
        
        
        if ~isempty(idx_band)
            features_depressed(i, j) = trapz(f_depressed(idx_band), pxx_depressed(idx_band));
        else
            features_depressed(i, j) = 0;   
        end
    end
end

 
disp('Extracted Features for Healthy Subjects:');
disp(features_healthy);

disp('Extracted Features for Depressed Subjects:');
disp(features_depressed);
% Visualize band powers for healthy and depressed groups
figure;
bar([mean(features_healthy); mean(features_depressed)]'); % Mean of features for each group
set(gca, 'XTickLabel', {'Delta', 'Theta', 'Alpha', 'Beta'});
ylabel('Mean Band Power');
title('Mean Band Power for Healthy vs. Depressed Subjects');
legend('Healthy', 'Depressed');
grid on;
 
specific_channel = 1; 
bands = {'Delta', 'Theta', 'Alpha', 'Beta'};
band_freqs = [0.5 4; 4 8; 8 13; 13 30];  
signal_specific = signal_healthy(:, specific_channel);  
sampling_rate = 256;  


[pxx_specific, f_specific] = pwelch(signal_specific, [], [], [], sampling_rate);
 
figure;

for j = 1:length(bands)
    subplot(2, 2, j);  
 
    plot(f_specific, 10*log10(pxx_specific));  
    hold on;

  
    xline(band_freqs(j, 1), 'r--', bands{j}, 'LabelHorizontalAlignment', 'left', 'Label', bands{j});
    xline(band_freqs(j, 2), 'r--', '', 'LabelHorizontalAlignment', 'right');

   
    title([bands{j} ' Band']);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    xlim([0 40]);   
    grid on;
    hold off;
end

sgtitle(['Frequency Decomposition of EEG Signal - Channel ', num2str(specific_channel)]);
 