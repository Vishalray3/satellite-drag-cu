function index_vec = starlink_noThrust(sat_ID)
if sat_ID == 3165
    index_vec = 450:1940;  % 2100:3680
elseif sat_ID == 3166
    index_vec  = 450:2000;
elseif sat_ID == 3167
    index_vec  = 930:2700; % 100:699  
elseif sat_ID == 3169
    index_vec  = 35:1700;  % 1:30 are sparse data
elseif sat_ID == 3174
    index_vec  = 500:2200;    
elseif sat_ID == 3178
    index_vec  = 720:2400;   
elseif sat_ID == 3181
    index_vec  = 520:2000;   
elseif sat_ID == 3182
    index_vec  = 670:2050;   
elseif sat_ID == 3189
    index_vec  = 530:1970;   
elseif sat_ID == 3401
    index_vec  = 300:1740;  
elseif sat_ID == 3415
    index_vec  = 650:2350;  
end
end