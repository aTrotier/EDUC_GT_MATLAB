function img = read_image_h5(filename)
if nargin < 1
    [file, PATHNAME] = uigetfile('*.h5');
    filename = fullfile(PATHNAME,file);
end

S=hdf5info(filename);

if exist(filename, 'file')
    dset = ismrmrd.Dataset(filename, 'dataset');
else
    error(['File ' filename ' does not exist.  Please generate it.'])
end

% hdr = ismrmrd.xml.deserialize(dset.readxml);

disp(dset.fid.identifier)

S=hdf5info(filename);
attributes=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(1).Name;
dataset=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(2).Name;
header=S.GroupHierarchy(1).Groups(1).Groups(1).Datasets(3).Name;

img=hdf5read(filename,dataset);
img=squeeze(img);
try
    imagine(img);
end
end