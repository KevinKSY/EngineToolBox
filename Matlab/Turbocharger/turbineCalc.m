load('SAE_map_turb_ABB_TPL_B');

TRef = SAE_map_turb.T_ref + 273.15;
pOut = 1e5 + 300*9.81;
pIn = pOut * SAE_map_turb.PR;
RPM = SAE_map_turb.RPM*sqrt(TRef);

[~,hBTurb,~,~,~,~,~,~,~,~,~,~,KBTurb] = GetThdynCombGasZachV1(pIn,TRef,1,0.0683);
TATurbIs = TRef.*(SAE_map_turb.PR).^((1-KBTurb)./KBTurb);
[~,hATurbIs,~,~,~,~,~,~,~,~,~,~,~] = GetThdynCombGasZachV1(pOut,TATurbIs,1,0.0683);

Uc_D = pi*RPM/30./ sqrt(2*(hBTurb-hATurbIs));
RPMUnique = unique(RPM);
figure
for i = 1:length(RPMUnique);
    idx = RPM == RPMUnique(i);
    plot(Uc_D(idx),SAE_map_turb.eff(idx),'*');
    hold on
    pause
end;
