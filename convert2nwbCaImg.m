% Convert two-photon calcium imaging data to the NWB format in Matlab
%
% Run this script to convert two-photon calcium imaging data and associated
% intracellular physiology data generated at the University of Bristol
% (UoB) to the Neurodata Without Borders (NWB) file format. This script is
% explained in the accompanying Bristol GIN for Calcium Imaging Data
% tutorial available at
% https://dervinism.github.io/bristol-neuroscience-data-guide/tutorials/Bristol%20GIN%20for%20Calcium%20Imaging%20Data.html
%
% You can use this script to get an idea of how to convert your own
% optical imaging data to the NWB file format.


%% Record metadata
% Project (experiment) metadata
projectName = 'Intracellular Ca2+ dynamics during plateau potentials trigerred by Schaffer collateral stimulation';
experimenter = 'MU';
institution = 'University of Bristol';
publications = 'In preparation';
lab = 'Jack Mellor lab';
brainArea = 'Hippocampus CA1-2';
greenIndicator = 'Fluo5f';
redIndicator = 'Alexa594';

% Animal metadata
animalID = 'm1';
ageInDays = 100;
age = ['P' num2str(ageInDays) 'D']; % Convert to ISO8601 format: https://en.wikipedia.org/wiki/ISO_8601#Durations
strain = 'C57BL/6J';
sex = 'M';
species = 'Mus musculus';
weight = [];
description = '001'; % Animal testing order.

% Session metadata
startYear = 2020;
startMonth = 12;
startDay = 4;
startTime = datetime(startYear, startMonth, startDay);
year = num2str(startYear); year = year(3:4);
month = num2str(startMonth); if numel(month) == 1; month = ['0' month]; end
day = num2str(startDay); if numel(day) == 1; day = ['0' day]; end
sliceNumber = 2;
cellNumber = 1;
imagingRate = 1/21; % A single linescan duration is 1sec with 20sec period in between linescans
lineRate = 1000.; % A number of lines scanned in a second.
sessionID = [animalID '_' year month day '_s' num2str(sliceNumber) ...
  '_c' num2str(cellNumber)]; % mouse-id_time_slice-id_cell-id
sessionDescription = 'Single cell imaging in a slice combined with somatic current clamp recordings and the stimulation of Schaffer collaterals';
sessionNotes = ['201204 - Slice 2 Imaging 3 dendrite regions roughly in the SO SR and SLM regions   ' ...
  'Same line scan with same intensity of stim (2.3V) at different locations along the cell   ' ...
  'Ephys frames to match up with linescans   ' ...
  'Frames   ' ...
  '1-8 Bottom Den   ' ...
  '10-19 Middle Den   ' ...
  '22-28 Top Dendrite   ' ...
  'Missed a few of the imaging with the ephys so more Ephys traces than linescans   ' ...
  'by the end of the experiment the top or neuron started to bleb.'];

% Generate Matlab classes from NWB core schema files
generateCore;

% Assign NWB file fields
nwb = NwbFile( ...
  'session_description', sessionDescription, ...
  'identifier', sessionID, ...
  'session_start_time', startTime, ...
  'general_experimenter', experimenter, ... % optional
  'general_session_id', sessionID, ... % optional
  'general_institution', institution, ... % optional
  'general_related_publications', publications, ... % optional
  'general_notes', sessionNotes, ... % optional
  'general_lab', lab); % optional

% Create subject object
subject = types.core.Subject( ...
  'subject_id', animalID, ...
  'age', age, ...
  'description', description, ...
  'species', species, ...
  'sex', sex);
nwb.general_subject = subject;

%% Convert calcium imaging data
% Create optical channels
green_optical_channel = types.core.OpticalChannel( ...
  'description', ['green channel corresponding to ' greenIndicator], ...
  'emission_lambda', 516.);

red_optical_channel = types.core.OpticalChannel( ...
  'description', ['red channel corresponding to ' redIndicator], ...
  'emission_lambda', 616.);

% Create the imaging device object
device = types.core.Device(...
  'description', 'Two-photon microscope', ...
  'manufacturer', 'Scientifica');
nwb.general_devices.set('2P_microscope', device);

% Create imaging plane objects
imaging_plane_name = 'green_imaging_plane';
green_imaging_plane = types.core.ImagingPlane( ...
  'optical_channel', green_optical_channel, ...
  'description', 'The plane for imaging calcium indicator Fluo5f.', ...
  'device', types.untyped.SoftLink(device), ...
  'excitation_lambda', 810., ...
  'imaging_rate', imagingRate, ...
  'indicator', 'Fluo5f', ...
  'location', brainArea);
nwb.general_optophysiology.set(imaging_plane_name, green_imaging_plane);

imaging_plane_name = 'red_imaging_plane';
red_imaging_plane = types.core.ImagingPlane( ...
  'optical_channel', red_optical_channel, ...
  'description', 'The plane for imaging calcium indicator Alexa594.', ...
  'device', types.untyped.SoftLink(device), ...
  'excitation_lambda', 810., ...
  'imaging_rate', imagingRate, ...
  'indicator', 'Alexa594', ...
  'location', brainArea);
nwb.general_optophysiology.set(imaging_plane_name, red_imaging_plane);

% Load data
dendrite = 'Bot';
data1 = load(['..\Analysed\' year month day '__s' num2str(sliceNumber) ...
  'd' num2str(cellNumber) '_004_ED__1 ' dendrite 'den_Analysed.mat']);
dendrite = 'Mid';
data2 = load(['..\Analysed\' year month day '__s' num2str(sliceNumber) ...
  'd' num2str(cellNumber) '_004_ED__1 ' dendrite 'den_Analysed.mat']);
dendrite = 'Top';
data3 = load(['..\Analysed\' year month day '__s' num2str(sliceNumber) ...
  'd' num2str(cellNumber) '_004_ED__1 ' dendrite 'den_Analysed.mat']);

% Append raw linescans so they have equal widths
data1.Analysed_data.Flur5_denoised = appendLinescans(data1.Analysed_data.Flur5_denoised);
data1.Analysed_data.Alexa_denoised = appendLinescans(data1.Analysed_data.Alexa_denoised);
data2.Analysed_data.Flur5_denoised = appendLinescans(data2.Analysed_data.Flur5_denoised);
data2.Analysed_data.Alexa_denoised = appendLinescans(data2.Analysed_data.Alexa_denoised);
data3.Analysed_data.Flur5_denoised = appendLinescans(data3.Analysed_data.Flur5_denoised);
data3.Analysed_data.Alexa_denoised = appendLinescans(data3.Analysed_data.Alexa_denoised);
nFrames = [size(data1.Analysed_data.Flur5_denoised,1) ...
  size(data2.Analysed_data.Flur5_denoised,1) size(data3.Analysed_data.Flur5_denoised,1)];

% Add optical physiology data: Fluorescence
input.indicator = greenIndicator;
input.imagingPlane = green_imaging_plane;
input.imagingRate = imagingRate;
input.lineRate = lineRate;
input.data = data1.Analysed_data.Flur5_denoised;
input.dendriteID = 'bottom';
nwb = setTwoPhotonSeries(nwb, input);

input.data = data2.Analysed_data.Flur5_denoised;
input.dendriteID = 'middle';
nwb = setTwoPhotonSeries(nwb, input);

input.data = data3.Analysed_data.Flur5_denoised;
input.dendriteID = 'top';
nwb = setTwoPhotonSeries(nwb, input);

input.indicator = redIndicator;
input.imagingPlane = red_imaging_plane;
input.data = data1.Analysed_data.Alexa_denoised;
input.dendriteID = 'bottom';
nwb = setTwoPhotonSeries(nwb, input);

input.data = data2.Analysed_data.Alexa_denoised;
input.dendriteID = 'middle';
nwb = setTwoPhotonSeries(nwb, input);

input.data = data3.Analysed_data.Alexa_denoised;
input.dendriteID = 'top';
nwb = setTwoPhotonSeries(nwb, input);
clear input

% Add optical physiology data: Delta fluorescence
input.indicator = redIndicator;
input.imagingPlane = red_imaging_plane;
input.imagingRate = imagingRate;
input.lineRate = lineRate;
input.data = data1.Analysed_data.Calcium_deltaF';
input.dendriteID = 'bottom';
nwb = setDeltaFSeries(nwb, input);

input.data = data2.Analysed_data.Calcium_deltaF';
input.dendriteID = 'middle';
nwb = setDeltaFSeries(nwb, input);

input.data = data3.Analysed_data.Calcium_deltaF';
input.dendriteID = 'top';
nwb = setDeltaFSeries(nwb, input);
clear input

% Add images
neuron_image = types.core.RGBImage( ...
  'data', data1.Analysed_data.Neuron, ...  % required: [height, width, colour]
  'description', 'RGB image of the full neuron.');

bottom_dend_image = types.core.GrayscaleImage( ...
  'data', data1.Analysed_data.ROI_img, ...  % required: [height, width]
  'description', 'Grayscale image of the bottom dendrite.');

middle_dend_image = types.core.GrayscaleImage( ...
  'data', data2.Analysed_data.ROI_img, ...  % required: [height, width]
  'description', 'Grayscale image of the middle dendrite.');

top_dend_image = types.core.GrayscaleImage( ...
  'data', data3.Analysed_data.ROI_img, ...  % required: [height, width]
  'description', 'Grayscale image of the top dendrite.');

image_collection = types.core.Images( ...
  'description', 'A collection of neuron and dendrite images.');
image_collection.image.set('neuron_image', neuron_image);
image_collection.image.set('dendrite1_image', bottom_dend_image);
image_collection.image.set('dendrite2_image', middle_dend_image);
image_collection.image.set('dendrite3_image', top_dend_image);

nwb.acquisition.set('ImageCollection', image_collection);


%% Convert intracellular electrophysiology data
% Create the recording device object
device = types.core.Device( ...
  'description', 'Amplifier for recording current clamp data.', ...
  'manufacturer', 'Molecular Devices');
nwb.general_devices.set('Amplifier_Multiclamp_700A', device);

electrode = types.core.IntracellularElectrode( ...
  'description', 'A patch clamp electrode', ...
  'location', 'Cell soma in CA1-2 of hippocampus', ...
  'slice', ['slice #' num2str(sliceNumber)], ...
  'device', types.untyped.SoftLink(device));
nwb.general_intracellular_ephys.set('icephys_electrode', electrode);

% Add current clamp data
input.ephysTime = data1.Analysed_data.Ephys_Time;
input.nSweeps = nFrames;
input.data = data1.Analysed_data.Ephys_data;
input.imagingRate = imagingRate;
input.electrode = electrode;
input.dendriteID = 'bottom';
nwb = setCClampSeries(nwb, input);

input.ephysTime = data2.Analysed_data.Ephys_Time;
input.data = data2.Analysed_data.Ephys_data;
input.dendriteID = 'middle';
nwb = setCClampSeries(nwb, input);

input.ephysTime = data3.Analysed_data.Ephys_Time;
input.data = data3.Analysed_data.Ephys_data;
input.dendriteID = 'top';
nwb = setCClampSeries(nwb, input);

%% Save the converted NWB file
nwbExport(nwb, [sessionID '.nwb']);



%% Local functions
function truncatedLinescans = truncateLinescans(linescans) %#ok<*DEFNU>
% truncatedLinescans = truncateLinescans(linescans)
%
% Function takes linescans of unequal widths and truncates them to have the
% smallest width and outputs the truncated linescans as a 3D matrix with
% the first dimension corresponding to the linescan number.

% Find the smallest width
nScans = numel(linescans);
length = size(linescans{1},1);
minWidth = size(linescans{1},2);
for iScan = 2:nScans
  minWidth = min([minWidth size(linescans{iScan},2)]);
end

% Truncate the linescans
truncatedLinescans = zeros(nScans, length, minWidth);
for iScan = 1:nScans
  diff = size(linescans{iScan},2) - minWidth;
  truncatedLinescans(iScan,:,:) = linescans{iScan}(:,1+floor(diff/2):end-ceil(diff/2));
end
end


function appendedLinescans = appendLinescans(linescans)
% appendedLinescans = appendLinescans(linescans)
%
% Function takes linescans of unequal widths and appends them to have the
% the same width and outputs the appended linescans as a 3D matrix with
% the first dimension corresponding to the linescan number.

% Find the largest width
nScans = numel(linescans);
length = size(linescans{1},1);
maxWidth = size(linescans{1},2);
for iScan = 2:nScans
  maxWidth = max([maxWidth size(linescans{iScan},2)]);
end

% Append the linescans
appendedLinescans = nan(nScans, length, maxWidth);
for iScan = 1:nScans
  appendedLinescans(iScan,:,1:size(linescans{iScan},2)) = linescans{iScan};
end
end


function nwb = setTwoPhotonSeries(nwb, input)
% nwb = setTwoPhotonSeries(nwb, input)
%
% Function creates, names, and adds a two-photon series data to a given NWB
% file.
%
% Input: nwb - the NWB file object.
%        input - a structure variable with the following fields:
%         indicator - a string variable with the calcium indicator name
%                     (e.g., 'Fluo5f' or 'Alexa594').
%         imagingPlane - the imaging plane object corresponding to the
%                        indicator.
%         imagingRate - a scalar with the rate at which individual full
%                       linescans are obtained.
%         lineRate - a scalar wit the frequency of individual lines within
%                    a single linescan.
%         data - an n-dimensional matrix containing two-photon series
%                linescan data. The first dimension corresponds to time
%                (that is, individual linescans). The second dimension
%                corresponds to individual lines along the vertical
%                denditic segment. The third dimension corresponds to the
%                dendritic width. NaNs are appended at the end of linescans
%                to make them equal in width.
%         dendriteID - a string variable with the dendritic ID (i.e.,
%                     'bottom', 'middle', and 'top').
%
% Output: nwb - the NWB file object containing the newly added two-photon
%               series data.

% Create image series
image_series = types.core.TwoPhotonSeries( ...
  'description', [input.indicator ' linescans of the ' input.dendriteID ' dendrite'], ...
  'imaging_plane', types.untyped.SoftLink(input.imagingPlane), ...
  'starting_time', 0.0, ...
  'starting_time_rate', input.imagingRate, ...
  'scan_line_rate', input.lineRate, ...
  'data', input.data, ...
  'data_unit', 'a.u.', ...
  'data_continuity', 'step', ...
  'comments', ['This two-photon series contains ' input.indicator ' linescans of the ', ...
  input.dendriteID ' (ROI) with the first dimension corresponding to time ', ...
  '(or to individual linescans). Each linescan is 1-sec in duration with ', ...
  '20-sec intervals between two linescans. The second dimension corresponds ', ...
  'to individual lines spanning the length of the dendrite in the ROI. ', ...
  'The third dimension corresponds to the width of the dendrite. ', ...
  'Some linescans may contain appended NaN values to make ', ...
  'widths of different linescans be equal within the same ROI.']);

% Name and add the two-photon series to the NWB file
if strcmpi(input.indicator, 'Fluo5f')
  opticalChannel = 'Green';
elseif strcmpi(input.indicator, 'Alexa594')
  opticalChannel = 'Red';
end
if strcmpi(input.dendriteID, 'bottom')
  dendriteID = '1';
elseif strcmpi(input.dendriteID, 'middle')
  dendriteID = '2';
elseif strcmpi(input.dendriteID, 'top')
  dendriteID = '3';
end
nwb.acquisition.set(['TwoPhotonSeries' opticalChannel dendriteID], image_series);
end


function nwb = setDeltaFSeries(nwb, input)
% nwb = setDeltaFSeries(nwb, input)
%
% Function creates, names, and adds a two-photon series delta fluorescence
% data to a given NWB file.
%
% Input: nwb - the NWB file object.
%        input - a structure variable with the following fields:
%         indicator - a string variable with the calcium indicator name
%                     (e.g., 'Fluo5f' or 'Alexa594').
%         imagingPlane - the imaging plane object corresponding to the
%                        indicator.
%         imagingRate - a scalar with the rate at which individual full
%                       linescans are obtained.
%         lineRate - a scalar wit the frequency of individual lines within
%                    a single linescan.
%         data - an 2D matrix containing two-photon series delta F data.
%                The first dimension corresponds to time (that is, individual
%                linescans). The second dimension corresponds to individual
%                lines along the vertical dendritic segment. The data is
%                averaged across the dendritic width.
%         dendriteID - a string variable with the dendritic ID (i.e.,
%                     'bottom', 'middle', and 'top').
%
% Output: nwb - the NWB file object containing the newly added two-photon
%               delta F series data.

% Create image series
image_series = types.core.TwoPhotonSeries( ...
  'description', ['Delta F data for the ' input.dendriteID, ...
  ' calculated based on ' input.indicator '.'], ...
  'imaging_plane', types.untyped.SoftLink(input.imagingPlane), ...
  'starting_time', 0.0, ...
  'starting_time_rate', input.imagingRate, ...
  'scan_line_rate', input.lineRate, ...
  'data', input.data, ...
  'data_unit', 'normalised', ...
  'data_continuity', 'step', ...
  'comments', ['This two-photon series contains delta F data calculated based on ', ...
  input.indicator ' for the ' input.dendriteID ' (ROI) ', ...
  'with the first dimension corresponding to time ', ...
  '(or to individual linescans). Each linescan is 1-sec in duration with ', ...
  '20-sec intervals between two linescans. The second dimension corresponds ', ...
  'to individual lines spanning the length of the dendrite in the ROI. ', ...
  'The data is averaged across the dendritic width.']);

% Name and add the two-photon delta F series to the NWB file
if strcmpi(input.dendriteID, 'bottom')
  dendriteID = '1';
elseif strcmpi(input.dendriteID, 'middle')
  dendriteID = '2';
elseif strcmpi(input.dendriteID, 'top')
  dendriteID = '3';
end
nwb.acquisition.set(['TwoPhotonDeltaFSeries' dendriteID], image_series);
end


function nwb = setCClampSeries(nwb, input)
% nwb = setCClampSeries(nwb, input)
%
% Function creates, names, and adds a current clamp series data to a given
% NWB file.
%
% Input: nwb - the NWB file object.
%        input - a structure variable with the following fields:
%         ephysTime - a time vector for a single recording sweep.
%         nSweeps - The total number of sweeps for every ROI.
%         imagingRate - a scalar with the rate at which individual full
%                       linescans are obtained.
%         data - a 2D matrix containing somatic current clamp recordings
%                corresponding to the initial part of the calcium
%                imaging period at dendritic ROI. The first dimension
%                corresponds to individual recording sweeps. The second
%                dimension corresponds to individual sweep data samples.
%         electrode - an electrode object.
%         dendriteID - a string variable with the dendritic ID (i.e.,
%                     'bottom', 'middle', and 'top').
%
% Output: nwb - the NWB file object containing the newly added current
%               clamp series data.

if strcmpi(input.dendriteID, 'bottom')
  nSweeps = input.nSweeps(1);
  sweeps = (1:nSweeps) - 1;
  dendriteID = '1';
elseif strcmpi(input.dendriteID, 'middle')
  nSweeps = input.nSweeps(2);
  sweeps = (input.nSweeps(1)+1:input.nSweeps(1)+nSweeps) - 1;
  dendriteID = '2';
else
  nSweeps = input.nSweeps(3);
  sweeps = (input.nSweeps(1)+input.nSweeps(2)+1:input.nSweeps(1)+input.nSweeps(2)+nSweeps) - 1;
  dendriteID = '3';
end

%timestamps = repmat((frames./input.imagingRate)', [1 numel(input.ephysTime)]) + ...
%  repmat((input.ephysTime./1000)', [nSweeps 1]); % sec
timestamps = input.ephysTime'./1000; % sec

current_clamp_series = types.core.CurrentClampSeries( ...
  'description', ['Somatic current clamp recording corresponding to ', ...
  'the initial part of the calcium imaging period ', ...
  'at ' input.dendriteID ' dendrite.'], ...
  'data', input.data(:,sweeps+1)', ...
  'data_continuity', 'step', ...
  'data_unit', 'millivolt', ...
  'starting_time', 0., ...
  'starting_time_rate', input.imagingRate, ...
  'electrode', types.untyped.SoftLink(input.electrode), ...
  'stimulus_description', 'N/A', ...
  'timestamps', timestamps, ...
  'comments', ['Somatic current clamp recording corresponding to ', ...
  'the initial part of the calcium imaging period ', ...
  'at ' input.dendriteID ' dendrite.', ...
  'The first dimension corresponds to individual ', ...
  'recording sweeps. The second dimension corresponds to ', ...
  'individual sweep data samples. The associated timestamps ', ...
  'variable provides timestamps for the second dimension. ', ...
  'This variable has to be combined with starting_time and ', ...
  'starting_time_rate variables to get absolute timestamps ', ...
  'for each data point']);
nwb.acquisition.set(['CurrentClampSeries' dendriteID], current_clamp_series);
end