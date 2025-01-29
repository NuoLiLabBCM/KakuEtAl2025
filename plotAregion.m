function plotAregion(region,colorNum)
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

p = patch(isosurface(xv.*sc, yv.*sc, zv.*sc, an ,1));
p.FaceAlpha = 0.1;
p.FaceColor = colorNum;
p.LineStyle = 'none'; 
set(gca, 'Color', [1 1 1], 'ZDir', 'Reverse');
% set(gca,'FontSize',18)
% view(51, 14); xlim([0 11.4]); ylim([0 13.2]); zlim([0 8]) % Allen