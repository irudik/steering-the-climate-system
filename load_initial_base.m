% Load initial guesses for base runs and all parameter combinations in the
% paper.

phi_base = .0091;
r_base = .055;
phi_h = phi_base*2;
phi_l = phi_base/2;
r_l = .014;

switch params.phi
    case phi_l
        switch params.r
            case r_l
                switch params.terminal_temp
                    case 2
                        terminal_search = [1339.2601];
                    case 2.5
                        terminal_search = [782.9722];
                    case 3
                        terminal_search = [244.6298];
                end
            case r_base
                switch params.terminal_temp
                    case 2
                        terminal_search = [9521.8895];
                    case 2.5
                        terminal_search = [5745.63];
                    case 3
                        terminal_search = [1812.621];
                end
        end
        
    case phi_base
        switch params.r
            case r_l
                switch params.terminal_temp
                    case 2
                        terminal_search = [675.1456];
                    case 2.5
                        terminal_search = [409.5148];
                    case 3
                        terminal_search = [133.9101];
                end
            case r_base
                switch params.terminal_temp
                    case 2
                        terminal_search = [4.031190900000000e+03];
                    case 2.5
                        terminal_search = [2584];
                    case 3
                        terminal_search = [878.4687];
                end
        end
    case phi_h
        switch params.r
            case r_l
                switch params.terminal_temp
                    case 2
                        terminal_search = [333.3439];
                    case 2.5
                        terminal_search = [202.5096];
                    case 3
                        terminal_search = [66.7691];
                end
            case r_base
                switch params.terminal_temp
                    case 2
                        terminal_search = [1604.944];
                    case 2.5
                        terminal_search = [1026.0613];
                    case 3
                        terminal_search = [351.9073];
                end
        end
end