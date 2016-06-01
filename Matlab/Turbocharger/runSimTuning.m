clear
eng.engLoad0 = 1;
EngineSystemParameters7X82;
mdl = 'TurbochargerFullMap';
open_system(mdl);
for ii = 1:11;
    eng.engLoad0 = eng_data.TC.Pe(ii)/eng_data.TC.Pe(2);
    EngineSystemParameters7X82;
    complete = 0;
    k = 0;
    while (complete == 0);
        k = k + 1;
        if k > 10
            break;
        end;
        try
            simOut = sim(mdl,'Solver','ode4','FixedStep','1e-3',...
                    'ZeroCross','on',...
                    'StopTime', '100', ... 
                    'ZeroCross','on', ...
                    'SaveTime','on',...
                    'SignalLogging','on','SignalLoggingName','logsout');
        catch err
            fprintf([err.identifier '\n']);
            close_system(mdl);
            open_system(mdl);
            continue
        end;
        complete = 1;
    end;
    answ = 'y';
    if k > 10
        propansw = 0;
        while propansw == 0;
            answ = input('Do you want to continue?(y/n) >','s');
            if answ == 'y' || answ == 'n'
                propansw = 1;
            end;
        end;
    end;
    if answ == 'n'
        break;
    end;
    logsout = simOut.get('logsout');
    omegaTC(ii) = logsout.get('omegaTC').Values.Data(end);


end;
        
    