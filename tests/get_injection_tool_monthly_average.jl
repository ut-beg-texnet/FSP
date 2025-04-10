using Dates
using CSV
using LinearAlgebra
using DataFrames  
using Statistics  


# Define the path to the CSV file containing well data
csv_name = "injection_tool_data.csv"

# Read the CSV file into a DataFrame structure
df = CSV.read(csv_name, DataFrame)

# Convert string dates to Date objects
df[!, "Date of Injection"] = Date.(df[!, "Date of Injection"])

# Extract unique well identifiers from the "API Number" column
# The [!, "column name"] syntax is used to access columns with spaces in their names
well_ids = unique(df[!, "API Number"])


# Process each well's injection data individually
for well_id in well_ids
    # Filter the DataFrame to only include rows for the current well
    well_data = df[df[!, "API Number"] .== well_id, :]
    
    # Extract the injection volume data and date of injection for this specific well
    injection_rate = well_data[!, "Volume Injected (BBLs)"]
    injection_date = well_data[!, "Date of Injection"]
    
    # Create a DataFrame with dates and injection rates
    well_df = DataFrame(
        "Date" => injection_date,
        "InjectionRate" => injection_rate
    )
    
    # Add columns for year and month to facilitate grouping
    well_df[!, :Year] = year.(well_df.Date)
    well_df[!, :Month] = month.(well_df.Date)
    
    # Group by year and month, then calculate average
    monthly_averages = combine(groupby(well_df, [:Year, :Month]), 
        :InjectionRate => mean => :AverageInjectionRate)
    
    # Add start and end date columns for each month
    monthly_averages[!, :StartDate] = Date.(monthly_averages.Year, monthly_averages.Month, 1)
    monthly_averages[!, :EndDate] = lastdayofmonth.(monthly_averages.StartDate)
    
    # Reorder columns for better readability
    monthly_averages = monthly_averages[:, [:StartDate, :EndDate, :Year, :Month, :AverageInjectionRate]]
    
    # Sort by date
    sort!(monthly_averages, [:Year, :Month])
    
    # Print the monthly averages for the current well
    println("Monthly average injection rates for well $well_id:")
    println(monthly_averages)
    println("------------------------------")
end














