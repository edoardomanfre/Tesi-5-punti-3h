using CSV, DataFrames

function print_discharge_csv(discharge, MM3Step)
    weekly = sum(MM3Step .* discharge[2, :, :, :], dims = 3)
    CSV.write(
        "Discharge_weekly.csv",
        DataFrame(transpose(weekly[:, :, 1])),
        delim = ";",
        header = false,
    )
end
