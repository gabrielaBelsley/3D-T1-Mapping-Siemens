function [diffAscconvCell] = diffAscconv(headerInfo1,headerInfo2,header1Name,header2Name)
    %diffAscconv difference between two ASCONV headers
    %e.g.: used to check that the only difference between GRE EPI Alpha and 2Alpha is limited to the FA and voltages 
    %where voltage 2alpha should be double of voltage for alpha 

    %   Input:
    %       headerInfo1 - header for Img1
    %       headerInfo2 - header for Img2
    %       header1Name - Name of header for Img1
    %       header2Name - Name of header for Img2
    %   Output:
    %       diffAscconvCell - cell containing the differences between Img 1
    %       and Img 2 header
    
    % Copyright, Gabriela Belsley (gabriela.belsley@stx.ox.ac.uk), 2021

    diffAscconvCell{1,1} = 'ASCCONV Field';
    diffAscconvCell{1,2} = header1Name;
    diffAscconvCell{1,3} = header2Name;
    cnt = 1;
    NFieldsHeader1 = length(headerInfo1);
    NFieldsHeader2 = length(headerInfo2);
    NCompare = min(NFieldsHeader1,NFieldsHeader2);
    for iascconv = 1:NCompare
       % headerInfo1{iascconv,1}, headerInfo2{iascconv,1}
        if isequal(headerInfo1{iascconv,1}, headerInfo2{iascconv,1})
            compareAscconv = isequal(headerInfo1{iascconv,2},headerInfo2{iascconv,2});
            
            if compareAscconv == 0
                cnt = cnt+1;
                
                % Name field that is different
                diffAscconvCell{cnt,1} = headerInfo1{iascconv,1};
                % Value for header 1
                diffAscconvCell{cnt,2} = headerInfo1{iascconv,2};
                % Value for header 2
                diffAscconvCell{cnt,3} = headerInfo2{iascconv,2};
                
                % Print to command Window
                fprintf('%s: ',headerInfo1{iascconv,1})
                if isnumeric(headerInfo1{iascconv,2})
                    fprintf('%0.4f, ',headerInfo1{iascconv,2})
                    fprintf('%0.4f \n',headerInfo2{iascconv,2})
                else
                    fprintf('%s, ',headerInfo1{iascconv,2})
                    fprintf('%s \n',headerInfo2{iascconv,2})
                end
                
            end
        end
    end
    
end

