# read price trajectory

function ReadPrice(path)
    f = open(path)
    Price = zeros(Float64, NStage, NStep)
    for iStage = 1:NStage
        #line = readline(f); items = split(line, " "); Price[iStage,1:NStep] = fill(parse(Float64, items[1]), NStep);
        line = readline(f)
        items = split(line, " ")
        for iStep = 1:NStep
            Price[iStage, iStep] = parse(Float64, items[iStep])
        end
    end
    return Price
end
