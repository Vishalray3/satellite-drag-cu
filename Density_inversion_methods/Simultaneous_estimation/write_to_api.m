function [] = write_to_api(disp_message, flag_disp)
if flag_disp == 1
url = 'https://lab.app.kayhan.io/api/demo-logger/write_log';
options = weboptions('HeaderFields',{'Content-Type' 'application/json'},'ContentType','json');
str_msg = strcat("{","'time'",": ", "'",string(datetime('now', 'TimeZone','Z', 'Format','yyyy-MM-dd''T''HH:mm:ss')),"'",", 'provider'",": ", "'Kayhan'", ", 'message'",": ",...
    disp_message,...
    ", 'wait'",": ", "'True'","}");
Sjson = struct('log', str_msg);
s = jsonencode(Sjson);
webwrite(url, s, options)
end