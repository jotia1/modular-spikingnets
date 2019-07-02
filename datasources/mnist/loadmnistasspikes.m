function [ spike_digits, digit_labels ] = loadmnistasspikes()
%% Assumes mnist_train.mat exists
load('mnist_train.mat', 'train_X', 'train_labels')

%% Filter only 0 and 8's
digit_0 = train_X(train_labels == 1, :);
digit_8 = train_X(train_labels == 9, :);
all_digits = [digit_0 ; digit_8];
digit_labels = [ones(size(digit_0, 1), 1); 9 * ones(size(digit_8, 1), 1) ];

%% Reshape
images = zeros(size(digit_labels, 1), 28, 28);
for i = 1 : size(all_digits, 1)
    images(i, :, :) = flipud(rot90(reshape(all_digits(i, :), [28, 28])));
end

%% Plot of sanity
do_plot = false;   %% WARNING plotting is excessive ... 
if do_plot
    figure;
    for i = 1:20
        subplot(4,5,i);
        imshow(squeeze(images(randi(size(images, 1), 1), :, :)));
    end
end
    
%% Convert to spike patterns
spike_digits = {};
for i = 1 : size(all_digits, 1)
    digit = squeeze(images(i, :, :));
    [spike_idxs, spike_times] = ind2sub([28, 28], find(flipud(digit) > 0.5));
    
    if do_plot
        figure
        subplot(1, 2, 1)
        plot(spike_times, spike_idxs, '.k');
        axis([0, 28, 0, 28]);
        subplot(1, 2, 2)
        imshow(digit);
        title('done')
    end
    
    spike_digits{end + 1} = [spike_idxs, spike_times];
    
end

end


