
function GT_reco_nonCart2(connection)
tic
disp("handle_connection was called.")

%connection.add_reader(uint32(gadgetron.Constants.RECON_DATA), @gadgetron.external.readers.read_recon_data_and_separated_density_weights);

next_acquisition = @connection.next;

acquisition = next_acquisition(); % Call input function to produce the next acquisition.

%% show traj and dcw
figure; 
subplot(2,1,1);plot(acquisition.bits.buffer.trajectory(1,:,1),...
										acquisition.bits.buffer.trajectory(2,:,1));
subplot(2,1,2);plot(acquisition.bits.buffer.trajectory(3,:,1));
toc
end

%% suppport functions
function reference = reference_header(recon_data)
    % We pick the first header from the header arrays - we need it to initialize the image meta data.    
    reference = structfun(@(arr) arr(:, 1)', recon_data.bits.buffer.headers, 'UniformOutput', false);
end
