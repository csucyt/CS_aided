% 该函数实现SNR和dB之间的转换
function result = dBChange(value,tag)

% value   : 输入的值
% tag     : 'dB2SNR' or 'SNR2dB'

if upper(tag(1)) == 'D'
    result = 10^(value/10);
elseif upper(tag(1)) == 'S'
    result = 10*log10(value);
else
    error('请输入db2SNR or SNR2db');
end

end