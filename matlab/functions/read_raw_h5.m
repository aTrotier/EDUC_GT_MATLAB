clear all

addpath('/home/ums/Downloads/NIfTI_20140122/')



% file='/data/mp2rage/Data/WHOBO_016_003/Original/WHOBO_016_003-3DT1_MP2RAGE08MM3_T1_orig.nii.gz';


% nifti_struct=load_nii(file);



%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% filename = '/data/dumpSiemens/Dernier_protocole/Semaine_du_080719/WHOBO16_MP2RAGE_090719_0_8mm.h5';
%filename = '/data/rawData_h5/WHOBO_19/WHOBO19_160719_MP2RAGE_0_8mm.h5'
 filename = '/data/dumpSiemens/ISMRMRD_UTE_Data_66056_5049812_5049817_23_20200221-145436.h5'
%filename = '/data/dumpSiemens/WHOBO_16/meas_MID00024_FID04374_iso0_8_angle7_CS2_85x2_85.h5'
S=hdf5info(filename);
if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr = ismrmrd.xml.deserialize(dset.readxml);

%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parametersclear all

addpath('/home/ums/Downloads/NIfTI_20140122/')



% file='/data/mp2rage/Data/WHOBO_016_003/Original/WHOBO_016_003-3DT1_MP2RAGE08MM3_T1_orig.nii.gz';


% nifti_struct=load_nii(file);



%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Loading an existing file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% filename = '/data/dumpSiemens/Dernier_protocole/Semaine_du_080719/WHOBO16_MP2RAGE_090719_0_8mm.h5';
%filename = '/data/rawData_h5/WHOBO_19/WHOBO19_160719_MP2RAGE_0_8mm.h5'
 filename = '/data/dumpSiemens/ISMRMRD_UTE_Data_66056_5049812_5049817_23_20200221-145436.h5'
%filename = '/data/dumpSiemens/WHOBO_16/meas_MID00024_FID04374_iso0_8_angle7_CS2_85x2_85.h5'
S=hdf5info(filename);
if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read some fields from the XML header %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to check if optional fields exists before trying to read them

hdr = ismrmrd.xml.deserialize(dset.readxml);

%% Encoding and reconstruction information
% Matrix size
enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
rec_Nz = hdr.encoding.reconSpace.matrixSize.z;

% Field of View
enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;

% Number of slices, coils, repetitions, contrasts etc.
% We have to wrap the following in a try/catch because a valid xml header may
% not have an entry for some of the parameters

try
  nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nSlices = 1;
end

try 
    nCoils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nCoils = 1;
end

try
    nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    nReps = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end

% TODO add the other possibilities

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time, 
% but this is much faster for data sets that fit into RAM.
D = dset.readAcquisition();


%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isCalib = D.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING')
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end
meas  = D.select(firstScan:D.getNumber);
clear D;

try
  nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
catch
    nSlices = 1;
end

try 
    nCoils = hdr.acquisitionSystemInformation.receiverChannels;
catch
    nCoils = 1;
end

try
    nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
catch
    nReps = 1;
end

try
    nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1 + 1;
catch
    nContrasts = 1;
end

% TODO add the other possibilities

%% Read all the data
% Reading can be done one acquisition (or chunk) at a time, 
% but this is much faster for data sets that fit into RAM.
D = dset.readAcquisition();


%% Ignore noise scans
% TODO add a pre-whitening example
% Find the first non-noise scan
% This is how to check if a flag is set in the acquisition header
isCalib = D.head.flagIsSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING')
isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
firstScan = find(isNoise==0,1,'first');
if firstScan > 1
    noise = D.select(1:firstScan-1);
else
    noise = [];
end
meas  = D.select(firstScan:D.getNumber);
clear D;


%% Reconstruct images
% Since the entire file is in memory we can use random access
% Loop over repetitions, contrasts, slices
reconImages = {};
nimages = 0;
for rep = 1:nReps
    for contrast = 1:nContrasts
        for slice = 1:nSlices
            % Initialize the K-space storage array
            K = zeros(enc_Nx, enc_Ny, enc_Nz, nCoils);
            % Select the appropriate measurements from the data
            acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
                        & (meas.head.idx.repetition==(rep-1)) ...
                        & (meas.head.idx.slice==(slice-1)));
            for p = 1:length(acqs)
                ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                K(:,ky,kz,:) = meas.data{acqs(p)};
            end
            % Reconstruct in x
            K = fftshift(ifft(fftshift(K,1),[],1),1);
            % Chop if needed
            if (enc_Nx == rec_Nx)
                im = K;
            else
                ind1 = floor((enc_Nx - rec_Nx)/2)+1;
                ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
                im = K(ind1:ind2,:,:,:);
            end
            % Reconstruct in y then z
            im = fftshift(ifft(fftshift(im,2),[],2),2);
            if size(im,3)>1
                im = fftshift(ifft(fftshift(im,3),[],3),3);
            end
            
            % Combine SOS across coils
            im = sqrt(sum(abs(im).^2,4));
            
            % Append
            nimages = nimages + 1;
            reconImages{nimages} = im;
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%


% nifti_struct.hdr.hist

hdr.encoding.encodedSpace.matrixSize.z
hdr.encoding.encodedSpace.fieldOfView_mm.z
hdr.encoding.reconSpace.matrixSize.z
hdr.encoding.reconSpace.fieldOfView_mm.z

hdr.encoding.reconSpace.fieldOfView_mm.x
hdr.encoding.reconSpace.fieldOfView_mm.y

meas.head.patient_table_position(:,1)

meas.head.position(:,1)

meas.head.read_dir(:,1)
meas.head.phase_dir(:,1)
meas.head.slice_dir(:,1)

rot= get_rot( meas.head );



param.Position=meas.head.position(:,1);
% quat = RawFile.image.slicePos(4:7,1);
a=quat(1);b=quat(2);c=quat(3);d=quat(4);
rot=[a*a+b*b-c*c-d*d,2*b*c-2*a*d,2*b*d+2*a*c;2*b*c+2*a*d,a*a+c*c-b*b-d*d,2*c*d-2*a*b;2*b*d-2*a*c,2*c*d+2*a*b,a*a+d*d-c*c-b*b];


% R11=meas.head.read_dir(1,1)
% R12=meas.head.read_dir(2,1)
% R13=meas.head.read_dir(3,1)
% 
% R21=meas.head.read_dir(1,2)
% R22=meas.head.read_dir(2,2)
% R23=meas.head.read_dir(3,2)
% 
% R31=meas.head.read_dir(1,3)
% R32=meas.head.read_dir(2,3)
% R33=meas.head.read_dir(3,3)

     a = 0.5  * sqrt(1+R11+R22+R33)   % (not stored)
     b = 0.25 * (R32-R23) / a      %  => quatern_b
     c = 0.25 * (R13-R31) / a      % => quatern_c
     d = 0.25 * (R21-R12) / a     %  => quatern_d
