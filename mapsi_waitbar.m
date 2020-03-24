function mapsi_waitbar(status, arg)

persistent numiter
persistent completed
persistent prevdone
persistent prevdone2
persistent warnings

switch status
    case 0 % setup
        warnings = [];
        
        fprintf('\n')
        
        tic
    case 1
        fprintf([arg '\n'])
    case 2
        start_bar(arg);
    case 3 % update
        if nargin == 1
            update_bar;
        else
            update_bar_comp(arg);
        end
    case 4 % finish
        if ~isempty(warnings)
           fprintf(['Warnings:\n', warnings, '\n']); 
        end
        
        tottime = ceil(toc);
        h = floor(tottime/3600);
        tottime = tottime - h*3600;
        m = floor(tottime/60);
        tottime = tottime - m*60;
        
        fprintf('Total time: %dh %dm %ds\n\n', h, m, tottime);
    case -1 % Warnings
        warnings = [warnings, arg, '\n'];
end

    function start_bar(tot)
        numiter = tot;
        completed = 0;
        prevdone2 = 0;
        prevdone = 0;
        
        fprintf(['[' repmat(' ',1,50) ']   0%%'])
    end

    function update_bar()
        completed = completed + 1;
        done = floor(completed/numiter*100);
        done2 = floor(completed/numiter*100 / 2);
        
        back = repmat('\b', 1, 6 + (50 - prevdone2));
        
        if done2 ~= prevdone2
            fprintf([back repmat('\x25A0',1,done2-prevdone2) repmat(' ',1,50-done2) '] %3d%%'], done);
            prevdone2 = done2;
            prevdone = done;
        elseif done ~= prevdone
            fprintf([repmat('\b',1,4) '%3d%%'], done)
            prevdone = done;
        end
        
        if completed == numiter
            fprintf('\nDone!\n\n')
        end
    end

    function update_bar_comp(completed)
        done = floor(completed/numiter*100);
        done2 = floor(completed/numiter*100 / 2);
        
        back = repmat('\b', 1, 6 + (50 - prevdone2));
        
        if done2 ~= prevdone2
            fprintf([back repmat('\x25A0',1,done2-prevdone2) repmat(' ',1,50-done2) '] %3d%%'], done);
            prevdone2 = done2;
            prevdone = done;
        elseif done ~= prevdone
            fprintf([repmat('\b',1,4) '%3d%%'], done)
            prevdone = done;
        end
        
        if completed == numiter
            fprintf('\nDone!\n\n')
        end
    end

end

