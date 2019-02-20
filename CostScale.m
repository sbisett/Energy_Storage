function [ TargetCost ] = CostScale(ScaleMethod,TargetSize,BaseSize,BaseCost,ScaleParameter)
%This code applies exponential scaling to the cost parameters
    if ScaleMethod == 1 %Exponential Method
        if isnan(ScaleParameter)
            ScaleParameter = 0.7;
        end
        TargetCost = BaseCost * (TargetSize/BaseSize)^ScaleParameter;
    elseif ScaleMethod == 2 %Linear Method
        TargetCost = BaseCost + (BaseCost/BaseSize)*(TargetSize-BaseSize);
    elseif ScaleMethod == 3 %Integer Method
        TargetCost = BaseCost * (TargetSize/BaseSize);
   % elseif ScaleMethod == 4 %Insulation
        
    else 
        error('ScaleMethod unrecognized')
    end

end

