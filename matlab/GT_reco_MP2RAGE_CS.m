
function GT_reco_MP2RAGE_CS(connection)
disp("handle_connection was called.")
addpath(genpath('functions')); % need a try

next_acquisition = @connection.next;

sRO = connection.header.encoding.reconSpace.matrixSize.x;
sPE = connection.header.encoding.reconSpace.matrixSize.y;
sSE = connection.header.encoding.reconSpace.matrixSize.z;
nContrast = connection.header.encoding.encodingLimits.contrast.maximum;
maxDim = nContrast + 1;

acquisition = next_acquisition(); % Call input function to produce the next acquisition.

nCh = size(acquisition.data.data,2);

kdata=zeros(sRO,nCh,sPE*sSE,maxDim); % initialize matrixwith zero


% find indices of each kspace line
row = acquisition.data.header.kspace_encode_step_1 + 1;
col = acquisition.data.header.kspace_encode_step_2 + 1;
TI_idx = acquisition.data.header.contrast + 1;

% buffer the data at the right place into kdata (faster than using
% bucket_to_buffer because of the sparse matrix)
kdata(:,:,sub2ind([sPE,sSE,maxDim],row,col,TI_idx))=acquisition.data.data;
kdata = permute(kdata,[1 3 2 4]);
kdata=reshape(kdata,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]

% Repeat this for reference data

% ref is deprecated you should use reference instead
kref=zeros(sRO,nCh,sPE*sSE,maxDim);
row = acquisition.reference.header.kspace_encode_step_1 + 1;
col = acquisition.reference.header.kspace_encode_step_2 + 1;

TI_idx = acquisition.reference.header.contrast + 1;

kref(:,:,sub2ind([sPE,sSE,maxDim],row,col,TI_idx))=acquisition.reference.data;
kref = permute(kref,[1 3 2 4]);
kref=reshape(kref,sRO,sPE,sSE,nCh,1,[]); %% buffering the echo according to bart convetion [RO,E1,E2,CHA,MAP,CON]

if  1 % plot mask
    mask = zeros(size(kref,1),size(kref,2),size(kref,3));
    mask(abs(kdata(:,:,:,1)) > 0) = 1;
    figure;imshow(squeeze(mask(1,:,:)),[]);
end

%% prepare parameter for reconstruction
% field associated with the sequence a_MP2RAGEPhase from the CRMSB, Bordeaux, France.

struct_MP2RAGE.calibSize = connection.header.userParameters.userParameterDouble(3).value;
struct_MP2RAGE.ETL = connection.header.userParameters.userParameterLong(6).value;
struct_MP2RAGE.TI1 = connection.header.userParameters.userParameterLong(2).value;
struct_MP2RAGE.TI2 = connection.header.userParameters.userParameterLong(3).value;
struct_MP2RAGE.alpha1 = connection.header.sequenceParameters.flipAngle_deg(1);
struct_MP2RAGE.alpha2 = connection.header.sequenceParameters.flipAngle_deg(1);
struct_MP2RAGE.MP2RAGE_TR = connection.header.userParameters.userParameterLong(4).value;
struct_MP2RAGE.TR = connection.header.sequenceParameters.TR;

%% Estimate coil sensitivity
run('/home/CODE/bart/startup.m'); % change this line for your bart installation
sensitivity_coil_map=bart(['caldir ' num2str(struct_MP2RAGE.calibSize)], kdata(:,:,:,:,:,2));
%
%sensitivity_coil_map=bart(['ecalib -m2 -r ' num2str(struct_MP2RAGE.calibSize)], kref(:,:,:,:,:,2));

%% Parallel + CS
struct.bg_mult=1;

disp('-------------------------------------------------');
disp('********** bart(pics...)  **********');

img_combined_CS = zeros(size(kdata,1),size(kdata,2),size(kdata,3),maxDim);
for i = 1:maxDim
    img_combined_CS(:,:,:,i) = bart('pics -g -d5 -e -i 60 -R W:7:0:0.01', kdata(:,:,:,:,:,i), sensitivity_coil_map);
end

[struct_T1map] = MP2RAGE_LookUpTable(img_combined_CS,struct_MP2RAGE); %% comb then MP2RAGE
disp('-----------RECO WITH correction background---------------');

multiFactor=struct.bg_mult*mean(mean(mean(abs(img_combined_CS(1:end,end-10:end,end-10:end,2)))));
MP2RAGE_CS=real((conj(img_combined_CS(:,:,:,1)).*img_combined_CS(:,:,:,2)-multiFactor)./(abs(img_combined_CS(:,:,:,1)).^2+abs(img_combined_CS(:,:,:,2)).^2+2*multiFactor));

%%

% img = gadgetron.lib.fft.cifftn(kdata,[1 2 3]); alternative using gadgetron fft
%
% img = fftshift(ifft(ifft(ifft(fftshift(kdata),[],1),[],2),[],3));
% img = squeeze(sqrt(sum(abs(img).^2,4)));
%
% figure;subplot(2,1,1);
% imshow(img(:,:,end/2,1),[]);
% subplot(2,1,2);
% imshow(img(:,:,end/2,2),[]);
% pause(3);
%
%
% img_to_send = img / max(img(:))*4096;

%% Prepare image to send
img_to_send = abs(img_combined_CS);
img_to_send(:,:,:,3) = MP2RAGE_CS;
img_to_send(:,:,:,4) = abs(struct_T1map.T1map);
%% send image

for i=1:size(img_to_send,4)
    image = gadgetron.types.Image.from_data(img_to_send(:,:,:,i), reference_header(acquisition));
    image.header.image_type = gadgetron.types.Image.MAGNITUDE;
    
    image.header.matrix_size(1) = size(image.data,1); % bad aloc with that
    image.header.matrix_size(2) = size(image.data,2);
    image.header.matrix_size(3) = size(image.data,3);
    image.header.channels=1;
    
    disp("Sending image to client.");
    connection.send(image);
end
end

%% suppport functions
function reference = reference_header(acquisition)
% We pick the first header from the header arrays - we need it to initialize the image meta data.
reference = structfun(@(arr) arr(:, 1)', acquisition.data.header, 'UniformOutput', false);
end
