%% Init
clc; clear all
init;
%% fig 2d (Recording coverage)
units_mtl=v_oralfacial_analysis.JawTuning * v_ephys.Unit * v_histology.ElectrodeCCFPositionElectrodePosition * v_ephys.UnitStat;
[colorVarMI, colorVarPh, kuiper_p, mi_perm, ccf_unitxA, ccf_unityA, ccf_unitzA, unit_amp, unit_fr] = units_mtl.fetchn('modulation_index', 'preferred_phase', 'kuiper_test', 'di_perm', 'ccf_x', 'ccf_y', 'ccf_z','unit_amp','avg_firing_rate');
% plot surface volumes 
plotBrain
hold on
% plot irn
plotAregion(136,'m');
% plot FN
plotAregion(661,'y');
% plot Hyp
plotAregion(773,'g');
% plot N. amb
plotAregion(939,'g');
% plot N. amb
plotAregion(143,'g');
% plot trigeminal
plotAregion(621,'c');
% vestibular nucleus
plotAregion([209 217],'b') % superior

% query and fetch tracks
mtl_track=v_ephys.ProbeInsertion * v_experiment.Session & 'rig="RRig-MTL"';
track_key=mtl_track.fetch();
% for i_track=round(9*length(track_key)/10):round(9.2*length(track_key)/10)
for i_track=1:length(track_key)
    sites=v_histology.ElectrodeCCFPositionElectrodePosition & track_key(i_track);
    [ccf_unitxT, ccf_unityT, ccf_unitzT] = sites.fetchn('ccf_x', 'ccf_y', 'ccf_z');
    hold on
    plot3(ccf_unitxT/1000,ccf_unitzT/1000,ccf_unityT/1000,'r','LineWidth', 0.01)
end

% sagittal view
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
% coronal view
% view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on

%% fig 2e (Firing rates)
% plot surface volumes 
plotBrain
hold on
% rate
MIthr=0;
hold on
crameri('lajolla',32)
caxis([0.2 50])  
axis on
%% sagittal view
temp = ccf_unitxA;
ccf_unitxA(:) = 5000;
scatter3(ccf_unitxA(colorVarMI>MIthr)/1000,ccf_unitzA(colorVarMI>MIthr)/1000,ccf_unityA(colorVarMI>MIthr)/1000,6,unit_fr(colorVarMI>MIthr),'filled')

% contour FN sag
plotAregionC_local(661,7.1,[],'k');
% plot Hyp sag
plotAregionC_local(773,6,[],'k');
% plot trigeminal sag
plotAregionC_local(621,7.3,[],'k');
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% vestibular nucleus
plotAregionC_local([209 217],7.34,[],'k') % lateral
plotAregionC_local(136,6.66,[],'k') % irn
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
ccf_unitxA = temp;
%% coronal view
temp = ccf_unitzA;
ccf_unitzA(:) = 9000;
scatter3(ccf_unitxA(colorVarMI>MIthr)/1000,ccf_unitzA(colorVarMI>MIthr)/1000,ccf_unityA(colorVarMI>MIthr)/1000,6,unit_fr(colorVarMI>MIthr),'filled')

% contour FN cor
plotAregionC_local(661,[],10.9,'k');
% plot Hyp cor
plotAregionC_local(773,[],12.5,'k');
% plot trigeminal cor
plotAregionC_local(621,[],10.2,'k');
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
plotAregionC_local([209 217],[],11,'k') % lateral
plotAregionC_local(136,[],11.9,'k') % irn
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
colorbar
ccf_unitzA=temp;
%% fig 2f (Breathing activity map)
units_mtl=v_oralfacial_analysis.BreathingTuning * v_ephys.Unit * v_histology.ElectrodeCCFPositionElectrodePosition;
[colorVarMI, colorVarPh, ccf_unitxA, ccf_unityA, ccf_unitzA] = units_mtl.fetchn('modulation_index', 'preferred_phase', 'ccf_x', 'ccf_y', 'ccf_z');
% plot brain
plotBrain
hold on
% masseter tracing
load('.\data\JunData\Masseter.mat')
Masset=table2array(roi_table{1,1}(:,3:5));
Masset(:,1)=5.4-Masset(:,1);
Masset(:,3)=5.7-Masset(:,3);
idx=randsample(size(Masset,1),size(Masset,1));
[count_m, edges_m, mid, loc] = histcn(Masset,10);
count_m = permute(count_m,[1 3 2]);
hold on
dX_m=(edges_m{3}(2)-edges_m{3}(1))/2;
dY_m=(edges_m{1}(2)-edges_m{1}(1))/2;
dZ_m=(edges_m{2}(2)-edges_m{2}(1))/2;
% genio
Genio=dlmread('.\data\JunData\024retro_genio3.txt',',',2,2);
Genio=Genio(:,1:3);
Genio(:,1)=5.4-Genio(:,1);
Genio(:,3)=5.7-Genio(:,3);
idx=randsample(size(Genio,1),size(Genio,1));
[count_g, edges_g, mid, loc] = histcn(Genio,10);
count_g = permute(count_g,[1 3 2]);
hold on
dX_g=(edges_g{3}(2)-edges_g{3}(1))/2;
dY_g=(edges_g{1}(2)-edges_g{1}(1))/2;
dZ_g=(edges_g{2}(2)-edges_g{2}(1))/2;
% jaw tuning
MIthr=0;
hold on

crameri('lajolla',32)
caxis([0.4 1]) % breath 
axis on
%% sagittal view
temp = ccf_unitxA;
ccf_unitxA(:) = 5000;
sizeVarMI = (colorVarMI*8)+5; %scale dot size according to modulation index (dotsize range 5-13)
[~,idx] = sort(colorVarMI);
scatter3(ccf_unitxA(idx)/1000, ccf_unitzA(idx)/1000, ccf_unityA(idx)/1000,sizeVarMI(idx),colorVarMI(idx),'filled')

% IRN
plotAregionC_local(136,6.66,[],'k');
% NTS
plotAregionC_local(651,6.71,[],'k');
% plot trigeminal sag
plotAregionC_local(621,7.3,[],'k');
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% contour FN sag
plotAregionC_local(661,7.1,[],'k');
% NA
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% plot Hyp sag
plotAregionC_local(773,6,[],'k');
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
ccf_unitxA = temp;
%% coronal view
temp = ccf_unitzA;
ccf_unitzA(:) = 9000;
sizeVarMI = (colorVarMI*8)+5; %scale dot size according to modulation index (dotsize range 5-13)
[~,idx] = sort(colorVarMI);

scatter3(ccf_unitxA(idx)/1000, ccf_unitzA(idx)/1000, ccf_unityA(idx)/1000,sizeVarMI(idx),colorVarMI(idx),'filled')

% IRN
plotAregionC_local(136,[],11.9,'k');
% NTS
plotAregionC_local(651,[],12,'k');
% plot trigeminal cor
plotAregionC_local(621,[],10.2,'k');
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% contour FN cor
plotAregionC_local(661,[],10.9,'k');
% % NA
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% plot Hyp cor
plotAregionC_local(773,[],12.5,'k');
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
colorbar
ccf_unitzA = temp;
%% fig 2g (licking activity map)
units_mtl=v_oralfacial_analysis.JawTuning * v_ephys.Unit * v_histology.ElectrodeCCFPositionElectrodePosition;
[colorVarMI, colorVarPh, kuiper_p, mi_perm, ccf_unitxA, ccf_unityA, ccf_unitzA] = units_mtl.fetchn('modulation_index', 'preferred_phase', 'kuiper_test', 'di_perm', 'ccf_x', 'ccf_y', 'ccf_z');
% plot brain
plotBrain
hold on
% masseter tracing
load('.\data\JunData\Masseter.mat')
Masset=table2array(roi_table{1,1}(:,3:5));
Masset(:,1)=5.4-Masset(:,1);
Masset(:,3)=5.7-Masset(:,3);
idx=randsample(size(Masset,1),size(Masset,1));
[count_m, edges_m, mid, loc] = histcn(Masset,10);
count_m = permute(count_m,[1 3 2]);
hold on
dX_m=(edges_m{3}(2)-edges_m{3}(1))/2;
dY_m=(edges_m{1}(2)-edges_m{1}(1))/2;
dZ_m=(edges_m{2}(2)-edges_m{2}(1))/2;
% genio
Genio=dlmread('.\data\JunData\024retro_genio3.txt',',',2,2);
Genio=Genio(:,1:3);
Genio(:,1)=5.4-Genio(:,1);
Genio(:,3)=5.7-Genio(:,3);
idx=randsample(size(Genio,1),size(Genio,1));
[count_g, edges_g, mid, loc] = histcn(Genio,10);
count_g = permute(count_g,[1 3 2]);
hold on
dX_g=(edges_g{3}(2)-edges_g{3}(1))/2;
dY_g=(edges_g{1}(2)-edges_g{1}(1))/2;
dZ_g=(edges_g{2}(2)-edges_g{2}(1))/2;
% jaw tuning
MIthr=0;
hold on
crameri('lajolla',32)
caxis([0.4 1]) % jaw 
axis on
%% sagittal view
temp = ccf_unitxA;
ccf_unitxA(:) = 5000;
sizeVarMI = (colorVarMI*8)+5; %scale dot size according to modulation index (dotsize range 5-13)
[~,idx] = sort(colorVarMI);
scatter3(ccf_unitxA(idx)/1000, ccf_unitzA(idx)/1000, ccf_unityA(idx)/1000,sizeVarMI(idx),colorVarMI(idx),'filled')

% IRN
plotAregionC_local(136,6.66,[],'k');
% NTS
plotAregionC_local(651,6.71,[],'k');
% plot trigeminal sag
plotAregionC_local(621,7.3,[],'k');
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% contour FN sag
plotAregionC_local(661,7.1,[],'k');
% NA
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% plot Hyp sag
plotAregionC_local(773,6,[],'k');
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
ccf_unitxA = temp;
%% coronal view
temp = ccf_unitzA;
ccf_unitzA(:) = 9000;
sizeVarMI = (colorVarMI*8)+5; %scale dot size according to modulation index (dotsize range 5-13)
[~,idx] = sort(colorVarMI);
scatter3(ccf_unitxA(idx)/1000, ccf_unitzA(idx)/1000, ccf_unityA(idx)/1000,sizeVarMI(idx),colorVarMI(idx),'filled')

% IRN
plotAregionC_local(136,[],11.9,'k');
% NTS
plotAregionC_local(651,[],12,'k');
% plot trigeminal cor
plotAregionC_local(621,[],10.2,'k');
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% contour FN cor
plotAregionC_local(661,[],10.9,'k');
% NA
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% plot Hyp cor
plotAregionC_local(773,[],12.5,'k');
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
colorbar
ccf_unitzA = temp;
%% fig 2h (swallow activity map)
load('.\data\figure_2\swallow.mat');
load('.\data\figure_2\pvals_logical_and_bootstrap_test_part_1.mat') 
load('.\data\figure_2\pvals_logical_and_bootstrap_test_part_2.mat')
modulation_index = modulation_indexNEW;
pvals = pvalsNEW;
ccf_x = double(ephys.ccf_x); % ccf x coord
ccf_y = double(ephys.ccf_y); % ccf y coord
ccf_z = double(ephys.ccf_z); %ccf z coord

%for bootstrapped test
mi = mean(modulation_index,2).';
num_events = active_trials;
pvals = max(pvals, [], 2).';

mi(pvals>0.01)=0; %change mod index for all non-sig neurons
mi(num_events<20)=0; %change mod index for all neurons less than 20 swallows
mi(mi<0.3) = 0;

colorVarMI = mi;

plotBrain
hold on
% masseter tracing
load('.\data\JunData\Masseter.mat')
Masset=table2array(roi_table{1,1}(:,3:5));
Masset(:,1)=5.4-Masset(:,1);
Masset(:,3)=5.7-Masset(:,3);
idx=randsample(size(Masset,1),size(Masset,1));
[count_m, edges_m, mid, loc] = histcn(Masset,10);
count_m = permute(count_m,[1 3 2]);
hold on
dX_m=(edges_m{3}(2)-edges_m{3}(1))/2;
dY_m=(edges_m{1}(2)-edges_m{1}(1))/2;
dZ_m=(edges_m{2}(2)-edges_m{2}(1))/2;
% genio
Genio=dlmread('.\data\JunData\024retro_genio3.txt',',',2,2);
Genio=Genio(:,1:3);
Genio(:,1)=5.4-Genio(:,1);
Genio(:,3)=5.7-Genio(:,3);
idx=randsample(size(Genio,1),size(Genio,1));
[count_g, edges_g, mid, loc] = histcn(Genio,10);
count_g = permute(count_g,[1 3 2]);
hold on
dX_g=(edges_g{3}(2)-edges_g{3}(1))/2;
dY_g=(edges_g{1}(2)-edges_g{1}(1))/2;
dZ_g=(edges_g{2}(2)-edges_g{2}(1))/2;


%% sagittal view
temp = ccf_x;
ccf_x(:) = 5000;
sizeVarMI = (colorVarMI*13)+5; %scale dot size according to modulation index (dotsize range 5-18)
[~,idx] = sort(colorVarMI);

scatter3(ccf_x(idx)/1000, ccf_z(idx)/1000, ccf_y(idx)/1000,sizeVarMI(idx),colorVarMI(idx), 'filled') 
crameri('lajolla',32)
caxis([0, 0.8])
axis on

% IRN
plotAregionC_local(136,6.66,[],'k');
% NTS
plotAregionC_local(651,6.71,[],'k');
% plot trigeminal sag
plotAregionC_local(621,7.3,[],'k');
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% contour FN sag
plotAregionC_local(661,7.1,[],'k');
% NA
plotAregionC_local(939,7.14,[],'k');
plotAregionC_local(143,7.04,[],'k');
% plot Hyp sag
plotAregionC_local(773,6,[],'k');
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
ccf_x = temp;
%% coronal view
temp = ccf_z;
ccf_z(:) = 9000;
sizeVarMI = (colorVarMI*13)+5; %scale dot size according to modulation index (dotsize range 5-18)
[~,idx] = sort(colorVarMI);

scatter3(ccf_x(idx)/1000, ccf_z(idx)/1000, ccf_y(idx)/1000,sizeVarMI(idx),colorVarMI(idx), 'filled')

crameri('lajolla',32)
caxis([0 0.8])
axis on

% IRN
plotAregionC_local(136,[],11.9,'k');
% NTS
plotAregionC_local(651,[],12,'k');
% plot trigeminal cor
plotAregionC_local(621,[],10.2,'k');
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% contour FN cor
plotAregionC_local(661,[],10.9,'k');
% NA
plotAregionC_local(939,[],11.7,'k');
plotAregionC_local(143,[],12.2,'k');
% plot Hyp cor
plotAregionC_local(773,[],12.5,'k');
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
colorbar
ccf_z = temp;
%%
function plotAregionC_local(region,xslice,yslice,color)
fn.AnnotatedBrain = 'F:\trackFinderData\Annotation_new_10_ds222_32bit.tif'; % 2017 v3

Anno = loadTifFast(fn.AnnotatedBrain);
if length(region)>1
    Anno(Anno~=region(1) & Anno~=region(2))=0;
else
    Anno(Anno~=region)=0;
end
sc = 0.02; %ccf
ds = 1;
an = Anno(1:ds:end, 1:ds:end, 1:ds:end);
an = permute(an,[3 2 1]);

xv = 1:ds:size(Anno, 2);
yv = 1:ds:size(Anno, 3);
zv = 1:ds:size(Anno, 1);

hold on

s=contourslice(xv.*sc, yv.*sc, zv.*sc, an,xslice,yslice,[]);
set(s,'EdgeColor',color,'LineWidth',0.5)
end
