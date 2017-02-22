function showmedida(name, val, ins, units)
    exponente10 = floor(log10(ins));
    if(ins < 2.5*10^exponente10)
        fprintf( [name ' = (%.1f \\pm %.1f) \\times 10^%d ' units '\n'], round(val/10^exponente10,1), round(ins/10^exponente10,1), exponente10);
    else
        fprintf( [name ' = (%.0f \\pm %.0f) \\times 10^%d ' units '\n'], round(val/10^exponente10), round(ins/10^exponente10), exponente10);
    end
end

