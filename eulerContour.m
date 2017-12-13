function eulerContour(fileBase, x0, x1, xh, y0, y1, yh)
    count = 0;
    anim = 1;
    while anim == 1
        %Each file would be labeled fileBase0, fileBase1, ..., so just check if
        %each file exists
        file = sprintf('%s%d', fileBase, count);
        if exist(file, 'file') == 2
            %File exists
            %Read in mass matrix
            mass = importdata(file);
            mass = flipud(mass);
            %Make contour
            X = linspace(x0, x1, xh);
            Y = linspace(y0, y1, yh);
            contour(X, Y, mass, 1);
            saveas(gcf, file, 'png')
            close
        else
            anim = 0;
        end
        count = count + 1;
    end
end