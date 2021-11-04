
function double = raw2double( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 16);
double = hex2num(hexstr);
end