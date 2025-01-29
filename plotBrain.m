function plotBrain
% plot surface volumes 
fn.AnnotatedBrain = '.\data\Annotation_new_10_ds222_32bit.tif'; % 2017 v3
Anno = loadTifFast(fn.AnnotatedBrain);
sc = 0.02; %ccf
ds = 1;
an = Anno(1:ds:end, 1:ds:end, 1:ds:end);
an = permute(an,[3 2 1]);

xv = 1:ds:size(Anno, 2);
yv = 1:ds:size(Anno, 3);
zv = 1:ds:size(Anno, 1);

% close all
figure; 
hold on
subplot(111); axis image;

p = patch(isosurface(xv.*sc, yv.*sc, zv.*sc, an ,1));
p.FaceAlpha = 0.05;
p.FaceColor = [0 0 0];
p.LineStyle = 'none';

set(gca, 'Color', [1 1 1], 'ZDir', 'Reverse');
set(gca,'FontSize',18)
set(gca,'TickDir','out')

view(180, 0); xlim([3 8.4]); ylim([9 13.2]); zlim([4 8])
axis off