% Figures for Chapter 3 Discussion 

addpath /Users/charlotteb/Documents/Chapter 3/Part 1/

load('Chapter3_results.mat');

% Eliminate scans which have different trend 
% < -80 for iFit (first three) 
% < -82 DOAS (first two) 

% shv = shv(4:end,:); 
% masaya1417 = masaya1417(3:end,:); 
% masaya1941 = masaya1941(3:end,:); 
% 
% masaya1417(:,5) = SO2D_1417(3:end,:); 
% masaya1941(:,5) = SO2D_1941(3:end,:); 

% Figure 14 
figure 
% ifit 
e = errorbar(shv(:,3), 1-shv(:,2), shv(:,4), shv(:,4), 'horizontal', 'or'); hold on; 
e.MarkerFaceColor = 'r';
e.MarkerEdgeColor = 'k';
% DOAS
% SCD
% masaya 1 
e = errorbar(masaya1417(:,3), 1-masaya1417(:,2), masaya1417(:,4), masaya1417(:,4), 'horizontal', 'ob'); 
e.MarkerFaceColor = 'b';
e.MarkerEdgeColor = 'k';
% masaya 2 
e = errorbar(masaya1941(:,3), 1-masaya1941(:,2), masaya1941(:,4), masaya1941(:,4), 'horizontal', 'oc'); 
e.MarkerFaceColor = 'c';
e.MarkerEdgeColor = 'k';
% % dSCD
% % masaya 1 
% e = errorbar(masaya1417(:,5), 1-masaya1417(:,2), masaya1417(:,4), masaya1417(:,4), 'horizontal', 'db'); 
% e.MarkerFaceColor = 'b';
% e.MarkerEdgeColor = 'k';
% % masaya 2 
% e = errorbar(masaya1941(:,5), 1-masaya1941(:,2), masaya1941(:,4), masaya1941(:,4), 'horizontal', 'dc'); 
% e.MarkerFaceColor = 'c';
% e.MarkerEdgeColor = 'k';
xlabel('SCD SO_2 (molec/cm^2)') 
ylabel('1 - minimum MSWC')
% legend('SHV (SCD, iFit)', 'Masaya 1417 (SCD, DOAS)', 'Masaya 1941 (SCD, DOAS)', 'Masaya 1417 (dSCD, DOAS)', 'Masaya 1941 (dSCD, DOAS)', 'Location','northwest');
% legend('SHV (iFit)', 'Masaya 14:17 UTC (DOAS)', 'Masaya 19:41 UTC (DOAS)', 'Location','northwest');
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
xlim([-0.2e18 3e18])
ylim([-0.02 1.02])
xline(shv(6,3), '-r');
xline(masaya1417(6,3), '-b');
xline(masaya1941(6,3), '-c');
legend('SHV (iFit)', 'Masaya 14:17 UTC (DOAS)', 'Masaya 19:41 UTC (DOAS)', 'Location','northwest');

% Clean up data 
% Convert to 1-MSWC
% ifit 
shv(:,2) = 1-shv(:,2); 
% DOAS 
masaya1417(:,2) = 1-masaya1417(:,2); 
masaya1941(:,2) = 1-masaya1941(:,2); 
% Remove low and high (background, saturation) 
% ifit 
for i = 1:length(shv)
    if (shv(i,2) < 0.005) || (shv(i,3) < 1.9e17)
       shv(i,2) = NaN;
    end 
end 
% DOAS
for i = 1:length(masaya1417)
    if (masaya1417(i,2) < 0.007) || (masaya1417(i,2) > 0.8) || (masaya1417(i,3) < 3e17)
       masaya1417(i,2) = NaN;
    end 
end 
for i = 1:length(masaya1941)
    if (masaya1941(i,2) < 0.008) || (masaya1941(i,3) < 2e17) 
       masaya1941(i,2) = NaN;
    end 
end 
    
% Figure 14 - clean 
figure 
% ifit 
e = errorbar(shv(:,3), shv(:,2), shv(:,4), shv(:,4), 'horizontal', 'or'); hold on; 
e.MarkerFaceColor = 'r';
e.MarkerEdgeColor = 'k';
% DOAS
% SCD
% masaya 1 
e = errorbar(masaya1417(:,3), masaya1417(:,2), masaya1417(:,4), masaya1417(:,4), 'horizontal', 'ob'); 
e.MarkerFaceColor = 'b';
e.MarkerEdgeColor = 'k';
% masaya 2 
e = errorbar(masaya1941(:,3), masaya1941(:,2), masaya1941(:,4), masaya1941(:,4), 'horizontal', 'oc'); 
e.MarkerFaceColor = 'c';
e.MarkerEdgeColor = 'k';
xlabel('SCD SO_2 (molec/cm^2)') 
ylabel('1 - minimum MSWC')
% legend('SHV (iFit)', 'Masaya 14:17 UTC (DOAS)', 'Masaya 19:41 UTC (DOAS)', 'Location','northwest');
Fig = gca; Fig.FontSize = 12; set(gcf,'color','w');
xlim([0 17e17])
ylim([0 0.8])
% xline(shv(6,3), '-r');
% xline(masaya1417(6,3), '-b');
% xline(masaya1941(6,3), '-c');
legend('SHV (iFit)', 'Masaya 14:17 UTC (DOAS)', 'Masaya 19:41 UTC (DOAS)', 'Location','northwest');

% Remove NaN 
shv(any(isnan(shv), 2), :) = [];  
masaya1417(any(isnan(masaya1417), 2), :) = [];  
masaya1941(any(isnan(masaya1941), 2), :) = [];  

% figure
[pshv, S] = polyfit(shv(:,3), shv(:,2), 1);
[fshv, deltashv] = polyval(pshv,shv(:,3), S); 
% plot(shv(:,3), shv(:,2),'or', shv(:,3),fshv,'-'); hold on 
plot(shv(:,3),fshv,':r'); hold on 
plot(shv(:,3),fshv+2*deltashv,'r-',shv(:,3),fshv-2*deltashv,'r-')

[pmasaya1417, S] = polyfit(masaya1417(:,3), masaya1417(:,2), 1);
[fmasaya1417, deltamasaya1417] = polyval(pmasaya1417,masaya1417(:,3), S); 
% plot(masaya1417(:,3), masaya1417(:,2),'ob', masaya1417(:,3),fmasaya1417,'-') 
plot(masaya1417(:,3),fmasaya1417,':b') 
plot(masaya1417(:,3),fmasaya1417+2*deltamasaya1417,'b-',masaya1417(:,3),fmasaya1417-2*deltamasaya1417,'b-')

[pmasaya1941, S] = polyfit(masaya1941(:,3), masaya1941(:,2), 1);
[fmasaya1941, deltamasaya1941] = polyval(pmasaya1941,masaya1941(:,3), S); 
% plot(masaya1941(:,3), masaya1941(:,2),'ob', masaya1941(:,3),fmasaya1941,'-') 
plot(masaya1941(:,3),fmasaya1941,':c') 
plot(masaya1941(:,3),fmasaya1941+2*deltamasaya1941,'c-',masaya1941(:,3),fmasaya1941-2*deltamasaya1941,'c-')

% Regression coefficient
% bVal is the slope or regression coefficient
% bVal = x\y
bValshv = masaya1417(:,3)\shv(:,2); 
bValmasaya1417 = masaya1417(:,3)\masaya1417(:,2); 
bValmasaya1941 = masaya1417(:,3)\masaya1941(:,2); 
% yCalc1 = b1 * x
yCalshv = bValshv*shv(:,3); 
yCalmasaya1417 = bValmasaya1417*masaya1417(:,3); 
yCalmasaya1941 = bValmasaya1941*masaya1941(:,3); 
% 
% figure
plot(shv(:,3), yCalshv, '-r'); hold on 
plot(masaya1417(:,3), yCalmasaya1417, '-b') 
plot(masaya1941(:,3), yCalmasaya1941, '-c') 
% Find x (SCD) by dividing  1-min(mswc)/bVal
bVals = [bValshv; bValmasaya1417; bValmasaya1941]; 
bValMean = mean(bVals); 














