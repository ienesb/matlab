function visualize_results(DD_map)
    % Visualize delay-Doppler map
    
    figure;
    imagesc(10*log10(fftshift(DD_map)));
    colorbar;
    title('Delay-Doppler Map (dB)');
    xlabel('Doppler Bins');
    ylabel('Delay Bins');
end
