%% SPC Data Plotting
% Uses the datasets assembled by dataImport.m and spcDataPull.m to plot
% figures from the SPC XP dataset.

%% Make Gravitational Correction
spc.d18Ograv = spc.d18O - 2*spc.d15N;
spc.d17Ograv = spc.d17O - 1*spc.d15N;
spc.d36Argrav = spc.d36Ar - 4*spc.d15N;
spc.d38Argrav = spc.d38Ar - 2*spc.d15N;
spc.dO2N2grav = spc.dO2N2 - 4*spc.d15N;
spc.dArN2grav = spc.dArN2 - 12*spc.d15N;

spc_replicates.d18Ograv = spc_replicates.d18O - 2*spc_replicates.d15N;
spc_replicates.d17Ograv = spc_replicates.d17O - 1*spc_replicates.d15N;
spc_replicates.d36Argrav = spc_replicates.d36Ar - 4*spc_replicates.d15N;
spc_replicates.d38Argrav = spc_replicates.d38Ar - 2*spc_replicates.d15N;
spc_replicates.dO2N2grav = spc_replicates.dO2N2 - 4*spc_replicates.d15N;
spc_replicates.dArN2grav = spc_replicates.dArN2 - 12*spc_replicates.d15N;

%% Calculate Pair Differences
pairDiffs=table;
pairDiffs.d15N = spc_replicates.d15N(:,1)-spc_replicates.d15N(:,2);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),1)-spc_replicates.d15N(isnan(pairDiffs.d15N),3);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),1)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),2)-spc_replicates.d15N(isnan(pairDiffs.d15N),3);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),2)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);
pairDiffs.d15N(isnan(pairDiffs.d15N),1) = spc_replicates.d15N(isnan(pairDiffs.d15N),3)-spc_replicates.d15N(isnan(pairDiffs.d15N),4);

pairDiffs.d18O = spc_replicates.d18O(:,1)-spc_replicates.d18O(:,2);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),1)-spc_replicates.d18O(isnan(pairDiffs.d18O),3);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),1)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),2)-spc_replicates.d18O(isnan(pairDiffs.d18O),3);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),2)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);
pairDiffs.d18O(isnan(pairDiffs.d18O),1) = spc_replicates.d18O(isnan(pairDiffs.d18O),3)-spc_replicates.d18O(isnan(pairDiffs.d18O),4);

pairDiffs.d17O = spc_replicates.d17O(:,1)-spc_replicates.d17O(:,2);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),1)-spc_replicates.d17O(isnan(pairDiffs.d17O),3);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),1)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),2)-spc_replicates.d17O(isnan(pairDiffs.d17O),3);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),2)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);
pairDiffs.d17O(isnan(pairDiffs.d17O),1) = spc_replicates.d17O(isnan(pairDiffs.d17O),3)-spc_replicates.d17O(isnan(pairDiffs.d17O),4);

pairDiffs.d36Ar = spc_replicates.d36Ar(:,1)-spc_replicates.d36Ar(:,2);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),1)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),1)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),2)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),2)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);
pairDiffs.d36Ar(isnan(pairDiffs.d36Ar),1) = spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),3)-spc_replicates.d36Ar(isnan(pairDiffs.d36Ar),4);

pairDiffs.d38Ar = spc_replicates.d38Ar(:,1)-spc_replicates.d38Ar(:,2);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),1)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),1)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),2)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),2)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);
pairDiffs.d38Ar(isnan(pairDiffs.d38Ar),1) = spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),3)-spc_replicates.d38Ar(isnan(pairDiffs.d38Ar),4);

pairDiffs.dO2N2 = spc_replicates.dO2N2(:,1)-spc_replicates.dO2N2(:,2);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),1)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),1)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),2)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),2)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);
pairDiffs.dO2N2(isnan(pairDiffs.dO2N2),1) = spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),3)-spc_replicates.dO2N2(isnan(pairDiffs.dO2N2),4);

pairDiffs.dArN2 = spc_replicates.dArN2(:,1)-spc_replicates.dArN2(:,2);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),1)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),3);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),1)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),2)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),3);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),2)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);
pairDiffs.dArN2(isnan(pairDiffs.dArN2),1) = spc_replicates.dArN2(isnan(pairDiffs.dArN2),3)-spc_replicates.dArN2(isnan(pairDiffs.dArN2),4);

pairDiffs.dO2N2grav = spc_replicates.dO2N2grav(:,1)-spc_replicates.dO2N2grav(:,2);
pairDiffs.dO2N2grav(isnan(pairDiffs.dO2N2grav),1) = spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),1)-spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),3);
pairDiffs.dO2N2grav(isnan(pairDiffs.dO2N2grav),1) = spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),1)-spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),4);
pairDiffs.dO2N2grav(isnan(pairDiffs.dO2N2grav),1) = spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),2)-spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),3);
pairDiffs.dO2N2grav(isnan(pairDiffs.dO2N2grav),1) = spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),2)-spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),4);
pairDiffs.dO2N2grav(isnan(pairDiffs.dO2N2grav),1) = spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),3)-spc_replicates.dO2N2grav(isnan(pairDiffs.dO2N2grav),4);

pairDiffs.dArN2grav = spc_replicates.dArN2grav(:,1)-spc_replicates.dArN2grav(:,2);
pairDiffs.dArN2grav(isnan(pairDiffs.dArN2grav),1) = spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),1)-spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),3);
pairDiffs.dArN2grav(isnan(pairDiffs.dArN2grav),1) = spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),1)-spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),4);
pairDiffs.dArN2grav(isnan(pairDiffs.dArN2grav),1) = spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),2)-spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),3);
pairDiffs.dArN2grav(isnan(pairDiffs.dArN2grav),1) = spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),2)-spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),4);
pairDiffs.dArN2grav(isnan(pairDiffs.dArN2grav),1) = spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),3)-spc_replicates.dArN2grav(isnan(pairDiffs.dArN2grav),4);

%% Pair Difference Plots
cluk; clc;

iBub = spc.bottomDepth < 700;
iBCTZ = spc.bottomDepth > 700 & spc.bottomDepth < 1140;
iClath = spc.bottomDepth > 1140;

figure; hold on;
P=gmregress(spc.dO2N2grav(iClath),spc.dArN2grav(iClath));
plot(spc.dO2N2grav,polyval(P,spc.dO2N2grav),'-k')
%plot(spc_replicates.dO2N2grav(iClath,:),spc_replicates.dArN2grav(iClath,:),'.','Color',lineCol(9));
plot(spc.dO2N2grav(iClath,:),spc.dArN2grav(iClath,:),'.','Color',lineCol(4));
xlabel('\deltaO_2/N_2 [per mil]'); ylabel('\deltaAr/N_2 [per mil]');
title('Ar/N_2 vs O_2/N_2 in SPC Clathrate Ice')

x = pairDiffs.dO2N2grav;
varsToPlot = [9 1:5]; numVars = length(varsToPlot);
depthsToPlot = [iBub iBCTZ any([iBub,iBCTZ],2) iClath]; numDepths = size(depthsToPlot,2);
variableNames = pairDiffs.Properties.VariableNames(varsToPlot);
massDiffs = [12 1 2 1 4 2]; cols = [lineCol(4); lineCol(1); lineCol(3); lineCol(3)*0.7; lineCol(5); lineCol(5)*0.7];
axCount = 1;

figure;
for ii = 1:numDepths
    for jj = 1:numVars
                
        y = pairDiffs{:,varsToPlot(jj)};
        
        [R,P] = corrcoef([x y],'Rows','Complete');
        r = R(1,2);
        p = P(1,2);
        
        
        subplot(numDepths,numVars,axCount); hold on;
        if p < 0.1
            [P,~,CI,pVal]=gmregress(x(depthsToPlot(:,ii)),y(depthsToPlot(:,ii)),0.001);
            plot(x,polyval(P,x),'-k');
            text(-8,0.75*min(get(gca,'YLim')),{['m = ' sprintf('%.2f',P(1)*1000) ' \pm ' sprintf('%.2f',(P(1)-CI(1,1))*1000) ' per meg/per mil'];['p = ' num2str(pVal)]});
        end
        
        plot(x(depthsToPlot(:,ii)),y(depthsToPlot(:,ii)),'.','Color',cols(jj,:))
        xlabel('\Delta\deltaO_2/N_2 [per mil]');
        ylabel(variableNames{jj});
        
        axCount=axCount+1;
    end
end
