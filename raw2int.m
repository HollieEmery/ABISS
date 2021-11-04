%% functions
function int = raw2int( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 8);
int = typecast(uint32(hex2dec(hexstr)),'int32');
end
