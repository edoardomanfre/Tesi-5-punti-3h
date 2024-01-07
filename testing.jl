function read_from_file_to_dict!(f, vars::Dict{Symbol, String})
    for line in eachline(f)
        if !isempty(strip(line)) && !startswith(strip(line),"#")
            var,val = strip.(split(line, "="))
            try
                vars[Symbol(strip(var))] = val
            catch end
        end
    end
    return vars
end

function read_to_dict(file, type::DataType)
    f = open(file)
    paramDict = Dict{Symbol, type}()
    paramDict =read_from_file_to_dict!(f, paramDict)

    return paramDict
end

function set_runCase(file = "runCase.in")
    
    runCaseDict = read_to_dict(file, String)
    case = caseData(;runCaseDict...)

    println(caseData)

    return caseData
end



paramDict = Dict{Symbol, String}()
paramDict = read_from_file_to_dict!("runCase.in", paramDict)

try 
    is(!true) 
catch e 
    println("error") 
end