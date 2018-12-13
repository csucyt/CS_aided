% �ú���ʵ��SNR��dB֮���ת��
function result = dBChange(value,tag)

% value   : �����ֵ
% tag     : 'dB2SNR' or 'SNR2dB'

if upper(tag(1)) == 'D'
    result = 10^(value/10);
elseif upper(tag(1)) == 'S'
    result = 10*log10(value);
else
    error('������db2SNR or SNR2db');
end

end