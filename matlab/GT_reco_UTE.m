
function GT_reco_UTE(connection)
tic
disp("handle_connection was called.")


next_acquisition = @connection.next;

sRO = connection.header.encoding.reconSpace.matrixSize.x;
sPE = connection.header.encoding.reconSpace.matrixSize.y;
sSE = connection.header.encoding.reconSpace.matrixSize.z;
nContrast = connection.header.encoding.encodingLimits.contrast.maximum;
maxDim = nContrast + 1;

acquisition = next_acquisition(); % Call input function to produce the next acquisition.

nCh = size(acquisition.data.data,2);

%% GENERAL PARAMETERS
addpath(genpath('/home/ums/Documents/mriSoft/gpuNUFFT-variousVer/gpuNUFFT'));
run('/home/ums/Documents/MATLAB/startup.m');
head = acquisition.data.header; %Just get header from first trajectory
matrice = acquisition.data.data;
matrice = permute(matrice,[1 3 2]);

ReadOffCenter  = connection.header.userParameters.userParameterDouble(1,5).value; % input for phase roll
PhaseOffCenter = connection.header.userParameters.userParameterDouble(1,6).value; % input for phase roll
acqOffSet      = connection.header.userParameters.userParameterDouble(1,7).value; % input for number of k=0 point (selfgating)
AcqMode        = connection.header.userParameters.userParameterDouble(1,8).value; % input : NONE ,REORDER, GODLEN ANGLE

useMultiCoil   = 1;
RecoMode       = 1; %1-Regridding 2-Regridding & sensibility map 3-CG-SENSE
SensitivityMap = 1; %1-2

N=2*size(matrice,1);
nSl=N;
[nRO,nPE,nCh]=size(matrice)



%% COMPUTE UTE traj
osf       = 1;  %oversampling
SGPoints  = 2*acqOffSet;
GA        = AcqMode;
RiseT     = 300;
BWp       = (1/(head.sample_time_us(1,1)*1e-6))/(osf*N); % Hz/px
traj      = ComputeUTEtraj(N, nPE , SGPoints, osf,BWp,RiseT, GA);

disp('traj')
disp(size(traj))

%% DCF
verbose = 1;
numIter = 10;
effMtx  = nRO;
DCF = sdc3_MAT(traj',numIter,effMtx,verbose,osf);

%% BART

%% bart reshape data

% 3D traj : [3,sizeR, proj]
traj = reshape(traj,sRO,[],3);
traj = permute(traj,[3,1,2]);

% kspace data : [1, sizeR, number of proj, number of channel]

matrice = reshape(matrice,1,sRO,nPE,nCh);

%% create DCF for each channel :
ss=size(matrice);
DCF_nCh=repmat(reshape(DCF,ss(1:3)),[1 1 1 nCh]);

% gridding

igrid = bart('nufft -i -c', traj*256,matrice(:,:,:,6));

igrid2 = bart('nufft -g -i -c', traj*256,matrice(:,:,:,6).*DCF_nCh);
igrid3 = bart('nufft -g -i -c', traj*256,matrice(:,:,:,6).*(DCF_nCh.^2));
igrid4 = bart('nufft -g -i -c', traj*256,matrice(:,:,:,6).*sqrt(DCF_nCh));
%% PHASE ROLL FOR OFFCENTER
for j=1:N
    for ii=1:N
        Offcenter(ii,j)= exp( (ii*2*pi*ReadOffCenter/(2*N) +j*2*pi*PhaseOffCenter/(2*N))*1i);
    end
end
for ii=1:N
    imgtemp(:,:,ii)=ifft2(fft2(img_comb(:,:,ii)).*Offcenter);
end

%% crop image
osf=2;
N=N/2;
image=abs(imgtemp(round((round(osf*N)-N)/2+1):round(N+(round(osf*N)-N)/2),round((round(osf*N)-N)/2+1):round(N+(round(osf*N)-N)/2),round((round(osf*N)-N)/2+1):round(N+(round(osf*N)-N)/2)));

%% IMAGE OUPUT
img_to_send=(abs(image)-min(min(min(abs(image)))))* 2^12/(max(max(max(abs(image))))-min(min(min(abs(image)))));

%% send image

for ii=1:size(img_to_send,4)
    image = gadgetron.types.Image.from_data(img_to_send(:,:,:,ii), reference_header(acquisition));
    image.header.image_type = gadgetron.types.Image.MAGNITUDE;
    
    image.header.matrix_size(1) = size(image.data,1); % bad aloc with that
    image.header.matrix_size(2) = size(image.data,2);
    image.header.matrix_size(3) = size(image.data,3);
    image.header.channels=1;
    
    disp("Sending image to client.");
    connection.send(image);
end
toc
end

%% suppport functions
function reference = reference_header(acquisition)
% We pick the first header from the header arrays - we need it to initialize the image meta data.
reference = structfun(@(arr) arr(:, 1)', acquisition.data.header, 'UniformOutput', false);
end
