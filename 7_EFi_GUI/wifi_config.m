function [PSDULength] = wifi_config(MCS)
    switch(MCS)
        case 0
            PSDULength = 100;
        case 1
            PSDULength = 200;
        case 2
            PSDULength = 300;
        case 3
            PSDULength = 400;
        case 4
            PSDULength = 600;
        case 5
            PSDULength = 800;
        case 6
            PSDULength = 900;
        case 7
            PSDULength = 1000;
    end    
end


