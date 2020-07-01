
function GT_pass_traj(connection)
tic
disp("handle_connection was called.")

next_acquisition = @connection.next;

flag = 1;
proj=0;
while flag
    proj=proj+1;
    try
        acquisition = next_acquisition(); % Call input function to produce the next acquisition.
    catch
        flag = 0;
    end
end
%% BART
% see webinard on how to use bart : https://github.com/mrirecon/bart-webinars


%% bart reshape data
% kspace data : [1, sRO, number of proj, number of channel]

sRO = size(matrice,1);
nCh = size(matrice,2);
NPro = size(matrice,3);

matrice2=permute(matrice,[1 3 2]);
matrice2=reshape(matrice2,1,sRO,[],nCh);
%% create traj
traj=compute_traj4bart(sRO, NPro, 0); % correspond to the uniform increment (360/NPro)

%% reco : 
igrid = bart('nufft -i', traj,matrice2); % iterative reconstruction 
igrid = bart('rss 8',igrid);
figure; subplot(1,3,1);
imshow(igrid,[]); title('Inverse');

%% Create traj with bart
traj_bart = bart(sprintf('traj -r -x %d -y %d',sRO,NPro));
igrid2 = bart('nufft -i', traj_bart,matrice2); % iterative reconstruction 
igrid2 = bart('rss 8',igrid2);
igrid2 = flipud(rot90(igrid2));
subplot(1,3,2);
imshow(igrid2,[]); title('Inverse');

%% Correct trajectory using ring

%GDring=bart('estdelay -R', traj_bart, matrice2);

%% Permute data : channel, readout, PE, SE
img_to_send{1}=permute(abs(igrid),[4, 1, 2, 3]);
img_to_send{2}=permute(abs(igrid2),[4, 1, 2, 3]);
%% send image

for ii=1:length(img_to_send)
    image = gadgetron.types.Image.from_data(img_to_send{ii}, acquisition.header);
    image.header.image_type = gadgetron.types.Image.MAGNITUDE;
    
    disp("Sending image to client.");
    connection.send(image);
end
toc
end

%% suppport functions

function K=compute_traj4bart(ADCres, Nspokes, SamplingType,TrajOffset)

% generate a radial trajectory with Nspokes lines.
% kloc_onesided=getpolar(Nspokes,ADCres);
% kloc_centered=kloc_onesided-ADCres/2-ADCres/2*1i-1-1i;

switch SamplingType
    case 0 % regular full spoke
        angleIncrement = pi / Nspokes;
    case 1 % golden angle full spoke
        angleIncrement = pi * (sqrt(5)-1)/2;
    case 2 % golden angle small version full spoke
        angleIncrement = pi * (3-sqrt(5))/2; 
    case 3 % regular half spoke
        angleIncrement = 2*pi / Nspokes;
    case 4 % golden angle half spoke
        angleIncrement = 2*pi * (sqrt(5)-1)/2;
    case 5% golden angle small version half spoke
        angleIncrement = 2*pi * (3-sqrt(5))/2; 
end

% between -.5 and .5
if(SamplingType<3)
   SpokeVector = linspace(-ADCres/2+1,ADCres/2,ADCres);
else
   SpokeVector = linspace(0,ADCres,ADCres);
end

% Compute the exact Fourier samples on the radial trajectory.

if(~exist('TrajOffset','var'))
    TrajOffset=0;
end

K = zeros([3,ADCres,Nspokes]);

 for s = 1:Nspokes
    cs = TrajOffset + s-1;
    K(1,:,s) = SpokeVector*cos(cs*angleIncrement);
    K(2,:,s) = SpokeVector*sin(cs*angleIncrement);

 end
end
