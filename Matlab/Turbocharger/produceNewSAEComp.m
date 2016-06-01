function [SAE_map_comp_new] = produceNewSAEComp(RPMNew,SAE_map_comp)
SAE_map_comp = redistSAEComp(SAE_map_comp);
mDot = SAE_map_comp.m_dot;
PR = SAE_map_comp.PR;
eff = SAE_map_comp.eff;
RPM = SAE_map_comp.RPM;

RPM_rep = unique(RPM);
m = length(RPM_rep);
ppMDot = cell(m,1);
ppEff = ppMDot;
for i = 1:m;
    idx = RPM == RPM_rep(i);
    mDotTemp = mDot(idx);
    PRTemp = PR(idx);
    effTemp = eff(idx);
    RPMTemp = RPM(idx);
    ppMDot = csaps(PRTemp,mDotTemp);
    ppEff = csaps(PRTemp,effTemp);
    PRGrid(i,:) = ppMDot.breaks;
    ppMDotGrid(i,:,:) = ppMDot.coefs;
    ppEffGrid(i,:,:) = ppEff.coefs;
    noOrd = ppMDot.order;
end;
% points interpolation
n = length(PRTemp);
for i = 1:n
    ppPR = csaps(RPM_rep,PRGrid(:,i));
    PRNew(:,i) = ppval(ppPR,RPMNew);
    if i < n
        for ii = 1:noOrd
            ppMDot = csaps(RPM_rep,ppMDotGrid(:,i,ii));
            ppEff = csaps(RPM_rep,ppEffGrid(:,i,ii));
            ppMDotGridNew(:,i,ii) = ppval(ppMDot,RPMNew);
            ppEffGridNew(:,i,ii) = ppval(ppEff,RPMNew);
        end;
    end;
end;
RPMNew = RPMNew*ones(1,n);
n = size(RPMNew);
for i = 1:n
    ppMDot.breaks = PRNew(i,:);
    ppMDot.coefs = permute(ppMDotGridNew(i,:,:),[2 3 1]);
    ppMDot.pieces = length(ppMDot.breaks) - 1;
    mDotNew(i,:) = ppval(ppMDot,PRNew(i,:));
    ppEff.breaks = PRNew(i,:);
    ppEff.coefs = permute(ppEffGridNew(i,:,:),[2 3 1]);
    ppEff.pieces = length(ppEff.breaks) - 1;
    effNew(i,:) = ppval(ppEff,PRNew(i,:));
end;
SAE_map_comp_new = SAE_map_comp;
SAE_map_comp_new.RPM = reshape(RPMNew',[],1);
SAE_map_comp_new.m_dot = reshape(mDotNew',[],1);
SAE_map_comp_new.eff = reshape(effNew',[],1);
SAE_map_comp_new.PR = reshape(PRNew',[],1);

function SAE_map_comp_new = redistSAEComp(SAE_map_comp)
% hObject    handle to compRedistPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SAE_map_comp_new = SAE_map_comp;
newNoPt = 10;
RPMComp_rep = num2cell(unique(SAE_map_comp.RPM));
newRPM = [];
newPR = [];
newMDot = [];
newEff = [];

for i = 1:length(RPMComp_rep);
    idx = find(SAE_map_comp.RPM == RPMComp_rep{i});
    tempRPM = SAE_map_comp.RPM(idx);
    tempPR = SAE_map_comp.PR(idx);
    tempMDot = SAE_map_comp.m_dot(idx);
    tempEff = SAE_map_comp.eff(idx);
    oldNoPt = length(tempRPM);
    rangePR = max(tempPR) - min(tempPR);
    rangeMDot = max(tempMDot) - min(tempMDot);
    rangeEff = max(tempEff) - min(tempEff);
    lArcPRMDot = sqrt(((tempPR(2:end)-tempPR(1:end-1))/rangePR).^2 + ...
        ((tempMDot(2:end)-tempMDot(1:end-1))/rangeMDot).^2);
    lArcPREff = sqrt(((tempPR(2:end)-tempPR(1:end-1))/rangePR).^2 + ...
        ((tempEff(2:end)-tempEff(1:end-1))/rangeEff).^2);
    totLArcPRMDot = sum(lArcPRMDot);
    newPRTemp = zeros(newNoPt,1);
    newMDotTemp = newPRTemp;
    newEffTemp = newPRTemp;
    newRPMTemp = newPRTemp;
    
    newPRTemp(1) = tempPR(1);
    newMDotTemp(1) = tempMDot(1);
    newEffTemp(1) = tempEff(1);
    newRPMTemp(1) = tempRPM(1);
    newPRTemp(newNoPt) = tempPR(oldNoPt);
    newMDotTemp(newNoPt) = tempMDot(oldNoPt);
    newEffTemp(newNoPt) = tempEff(oldNoPt);
    newRPMTemp(newNoPt) = tempRPM(oldNoPt);
    for iii=1:newNoPt-2
        k = 1;
        lArc = totLArcPRMDot/(newNoPt-1)*iii;
        while(lArc>sum(lArcPRMDot(1:k)))
            k = k + 1;
        end;
        newPRTemp(iii+1) = tempPR(k) + (tempPR(k+1) - tempPR(k))* ...
            (lArc - sum(lArcPRMDot(1:k-1)))/lArcPRMDot(k);
    end;
    pp1 = csaps(tempPR,tempMDot);
    pp2 = csaps(tempPR,tempEff);
    newMDotTemp = ppval(pp1,newPRTemp);
    newEffTemp = ppval(pp2,newPRTemp);
    newRPMTemp = tempRPM(1)*ones(newNoPt,1);
    newPR = [newPR; newPRTemp];
    newMDot = [newMDot; newMDotTemp];
    newEff = [newEff; newEffTemp];
    newRPM = [newRPM;newRPMTemp];
end;
SAE_map_comp_new.m_dot = newMDot;
SAE_map_comp_new.RPM = newRPM;
SAE_map_comp_new.PR = newPR;
SAE_map_comp_new.eff = newEff;