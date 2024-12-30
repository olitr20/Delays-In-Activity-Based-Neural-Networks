function c_palette = customColourPalette(colours)

    rgb = zeros(3,3);
    for i = 1:length(colours)
        r = hex2dec(colours{i}(1:2));
        g = hex2dec(colours{i}(3:4));
        b = hex2dec(colours{i}(5:6));
        rgb(i,:) = [r, g, b] / 255;
    end

    c_palette = [linspace(rgb(1,1), rgb(2,1), 256/2)', ...
                  linspace(rgb(1,2), rgb(2,2), 256/2)', ...
                  linspace(rgb(1,3), rgb(2,3), 256/2)'; ...
                  linspace(rgb(2,1), rgb(3,1), 256/2)', ...
                  linspace(rgb(2,2), rgb(3,2), 256/2)', ...
                  linspace(rgb(2,3), rgb(3,3), 256/2)'];
end