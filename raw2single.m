
function single = raw2single( raw )
hexstr = reshape(transpose(flip( dec2hex(raw, 2) )), 1, 8);
single = typecast(uint32(hex2dec(hexstr)),'single');
end
