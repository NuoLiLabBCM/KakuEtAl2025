%% Init
clc; clear all
cd C:\Users\NuoLiLabTower2\Documents\GitHub\map-ephys\MATLAB\pipeline
init;
%% fig 3a
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
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,7.26,[],[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,6.6,[],[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)
% plot irn
plotAregion(136,'m');
% plot FN
plotAregion(661,'y');
% plot Hyp
plotAregion(773,'g');
% plot trigeminal
plotAregion(621,'c');
% plot NA
plotAregion(939,'g');
plotAregion(143,'g');
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on
%% coronal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,[],10.9,[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,[],12,[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)
% plot irn
plotAregion(136,'m');
% plot FN
plotAregion(661,'y');
% plot Hyp
plotAregion(773,'g');
% plot trigeminal
plotAregion(621,'c');
% plot NA
plotAregion(939,'g');
plotAregion(143,'g');
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on

%% figure 3b
units_mtl=v_oralfacial_analysis.JawTuning * v_ephys.Unit * v_histology.ElectrodeCCFPositionElectrodePosition;
[colorVarMI, colorVarPh, kuiper_p, mi_perm, ccf_unitxA, ccf_unityA, ccf_unitzA] = units_mtl.fetchn('modulation_index', 'preferred_phase', 'kuiper_test', 'di_perm', 'ccf_x', 'ccf_y', 'ccf_z');

plotBrain
hold on

%plot jaw tuning
MIthr=0;
crameri('lajolla',32)
caxis([0.4 1]) % jaw 
axis on

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
temp = ccf_unitxA;
ccf_unitxA(:) = 5000;
sizeVarMI = (colorVarMI*8)+5; %scale dot size according to modulation index (dotsize range 5-13)
[~,idx] = sort(colorVarMI);
scatter3(ccf_unitxA(idx)/1000, ccf_unitzA(idx)/1000, ccf_unityA(idx)/1000,sizeVarMI(idx),colorVarMI(idx),'filled')


s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,7.26,[],[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,6.6,[],[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

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

s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,[],10.9,[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,[],12,[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

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

%% fig 3d (ALM axons)
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



load ./data/figure_3/analysis_anatomy_latALM_proj_medulla450
threshold = 0.6; %cut-off threshold for plotting projections
ds =1; %downsampling factor
sc = 0.02; %scaling factor

an = im_average(1:ds:end, 1:ds:end, 1:ds:end);
an = permute(an,[3 2 1]);

xv = 1:ds:size(im_average, 2);
yv = 1:ds:size(im_average, 3);
zv = 1:ds:size(im_average, 1);

p = patch(isosurface(xv(285:end).*sc, yv.*sc, zv.*sc, an(:,285:end,:) , threshold)); %only left hemisphere
p.FaceAlpha = 0.1;
p.FaceColor = 'g';
p.LineStyle = 'none'; 
set(gca, 'Color', [1 1 1], 'ZDir', 'Reverse');

%% sagittal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,7.26,[],[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,6.6,[],[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on
%% coronal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,[],10.9,[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,[],12,[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on

%% fig 3f (SC axons)
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

load ./data/figure_3/analysis_anatomy_latSC_proj_medulla450

threshold = 0.6; %cut-off threshold for plotting projections
ds =1; %downsampling factor
sc = 0.02; %scaling factor

an = im_average(1:ds:end, 1:ds:end, 1:ds:end);
an = permute(an,[3 2 1]);

xv = 1:ds:size(im_average, 2);
yv = 1:ds:size(im_average, 3);
zv = 1:ds:size(im_average, 1);

p = patch(isosurface(xv(285:end).*sc, yv.*sc, zv.*sc, an(:,285:end,:) , threshold)); %only left hemisphere
p.FaceAlpha = 0.1;
p.FaceColor = 'g';
p.LineStyle = 'none'; 
set(gca, 'Color', [1 1 1], 'ZDir', 'Reverse');

%% sagittal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,7.26,[],[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,6.6,[],[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on
%% coronal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,[],10.9,[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,[],12,[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
axis on

%% fig 3h
%fiber locations aligned (phox2b x ai32 only)
% x_coords = [6.46	6.45	6.53	6.29	6.96	4.85	7	4.4	6.8	4.66	6.42	6.7	4.8	6.35	4.6	5.02	4.6	5.02	6.72	4.46	6.67	4.78	6.64	4.77	6.22];
% y_coords = [11.98	11.39	11.51	11.41	11.7	11.7	11.07	10.97	11.15	10.86	12.26	11.71	11.62	12.26	12.03	11.85	11.46	11.51	11.4	11.4	12.15	12.07	12.02	11.95	12.31];
% z_coords = [6	6	5.3	5.41	5.12	5.27	6.05	6.11	5.32	5.84	5.38	5.74	6.1	5.8	5.93	5.62	5.8	5.1	5.66	5.2	5.05	5.1	5.87	5.1	5.8];
% colors = repmat(["b"], 1,size(x_coords,2));
% rhythmic_targets = [5,6,11,12,13,14,15,16,17,18,19,20,21,22,24,25];
% for r=rhythmic_targets
%     colors(r)='r';
% end

% fiber locations aligned (virus injection only)
% animal_ids = {'HK104_L'	'HK104_R'	'HK105_L'	'HK106_L'	'HK106_R'	'HK107_L'	'HK107_R'	'HK108_L'	'HK108_R'	'HK109_L'	'HK109_R'	'HK111_L'	'HK111_R'	'HK112_L'	'HK112_R'	'HK113_L'	'HK113_R'	'HK114_L'	'HK114_R'	'HK115_L'	'HK115_R'}
% x_coords = [6.27	4.73	6.25	6.2	4.52	6.47	4.7	5.14	5.05	7.27	4.72	7.23	4.18	6.24	3.73	6.99	3.6	5.4	4.05	6.96	6.42];
% y_coords = [11.69	12.14	11.02	11.12	11.68	11.19	12.21	10.72	11.45	9.96	12.4	11.5	11.1	11.58	10.48	10.02	10.91	12.85	11.88	10.26	11.71];
% z_coords = [5.59	6.07	5.46	5.8	5.77	5.67	5.88	5.92	5.14	5.36	5.62	5.66	4.85	5.82	4.63	5.44	5.2	5.71	5.95	5.58	5.36];
% colors = repmat(["b"], 1,size(x_coords,2));
% rhythmic_targets = [1,4,5,6,7,11,12,19];
% spc_targets = [10,13,15,17];
% for r=rhythmic_targets
%     colors(r)='r';
% end
% 
% for s=spc_targets
%     colors(s)='g';
% end

%fiber locations aligned (all animals)
animal_ids= {'HK7825'	'HK7827'	'HK7937_L'	'HK7938_L'	'HK8213_L'	'HK8213_R'	'HK7788_L'	'HK7788_R'	'HK7935_L'	'HK7935_R'	'HK8112'	'HK8120_L'	'HK8120_R'	'HK8280_L'	'HK8280_R'	'HK8281_L'	'HK8281_R'	'HK8447_R'	'HK8448_L'	'HK8448_R'	'HK8516_L'	'HK8516_R'	'HK8366_L'	'HK8366_R'	'HK8413_L'	'HK104_L'	'HK104_R'	'HK105_L'	'HK106_L'	'HK106_R'	'HK107_L'	'HK107_R'	'HK108_L'	'HK108_R'	'HK109_L'	'HK109_R'	'HK111_L'	'HK111_R'	'HK112_L'	'HK112_R'	'HK113_L'	'HK113_R'	'HK114_L'	'HK114_R'	'HK115_L'	'HK115_R'};
x_coords = [6.46	6.45	6.53	6.29	6.96	4.85	7	4.4	6.8	4.66	6.42	6.7	4.8	6.35	4.6	5.02	4.6	5.02	6.72	4.46	6.67	4.78	6.64	4.77	6.22	6.27	4.73	6.25	6.2	4.52	6.47	4.7	5.14	5.05	7.27	4.72	7.23	4.18	6.24	3.73	6.99	3.6	5.4	4.05	6.96	6.42];
y_coords = [11.98	11.39	11.51	11.41	11.7	11.7	11.07	10.97	11.15	10.86	12.26	11.71	11.62	12.26	12.03	11.85	11.46	11.51	11.4	11.4	12.15	12.07	12.02	11.95	12.31	11.69	12.14	11.02	11.12	11.68	11.19	12.21	10.72	11.45	9.96	12.4	11.5	11.1	11.58	10.48	10.02	10.91	12.85	11.88	10.26	11.71];
z_coords = [6	6	5.3	5.41	5.12	5.27	6.05	6.11	5.32	5.84	5.38	5.74	6.1	5.8	5.93	5.62	5.8	5.1	5.66	5.2	5.05	5.1	5.87	5.1	5.8	5.59	6.07	5.46	5.8	5.77	5.67	5.88	5.92	5.14	5.36	5.62	5.66	4.85	5.82	4.63	5.44	5.2	5.71	5.95	5.58	5.36];
colors = repmat(["b"], 1,size(x_coords,2));
rhythmic_targets = [5,6,11,12,13,14,15,16,17,18,19,20,21,22,24,25,26,29,30,31,32, 36,37,44];
spc_targets = [35,38,40,42];
for r=rhythmic_targets
    colors(r)='r';
end
for s=spc_targets
    colors(s)='g';
end

plotBrain
hold on
axis on

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
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,7.26,[],[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,6.6,[],[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

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
for i=1:length(x_coords)
    s(i) = scatter3(x_coords(i), y_coords(i),z_coords(i), 40,colors(i), "o","filled", 'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end
view(90, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom

%% coronal view
s=contourslice(edges_m{3}+dX_m, edges_m{1}+dY_m, edges_m{2}+dZ_m, count_m,[],10.9,[],[2 2]);
set(s,'EdgeColor','b','LineWidth',2)
s=contourslice(edges_g{3}+dX_g, edges_g{1}+dY_g, edges_g{2}+dZ_g, count_g,[],12,[],[6 6]);
set(s,'EdgeColor',[0.3 0.3 0.3],'LineWidth',2)

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
for i=1:length(x_coords)
    s(i) = scatter3(x_coords(i), 13.2,z_coords(i), 50,colors(i), "o","filled", 'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
end
view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([3 8]) % zoom
%%
function plotAregionC_local(region,xslice,yslice,color)
fn.AnnotatedBrain = '.\data\Annotation_new_10_ds222_32bit.tif'; % 2017 v3

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
set(s,'EdgeColor',color,'LineWidth',1)
end