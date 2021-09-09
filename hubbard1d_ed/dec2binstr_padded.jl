function dec2binstr_padded(i, NBitMax)
    bitStr = string(i,base=2)
    Npad = NBitMax - length(bitStr)
    return "0"^(Npad) * bitStr
end
# Convert from binstr to decimal
# parse(Int64, "1100", base=2)