function pow = coex_power( input )
    pow = sum(input.*conj(input))./length(input);
end
