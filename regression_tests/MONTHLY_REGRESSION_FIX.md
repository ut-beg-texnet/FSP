# Monthly Regression Test Fix

## Problem Summary

The `fsp_hydrology_monthly_regression.ipynb` notebook was showing differences between Julia and MATLAB results, particularly near the injection well. The differences appeared in the "Difference" heatmaps as red/blue spots close to the well location.

## Root Causes Identified

### 1. **Inconsistent Function Signatures**
The monthly notebook was using a different function signature than the working annual regression notebook:

**Before (Monthly - INCORRECT):**
```julia
function pressureScenario_Rall(
    bpds::Vector{Float64},
    days::Vector{Float64},
    r_meters::AbstractVector{Float64},
    STRho::Tuple{Float64, Float64, Float64},  # ❌ Tuple parameter
    evaluation_days_from_start::Union{Float64, Nothing} = nothing
)
```

**After (Monthly - CORRECTED):**
```julia
function pressureScenario_Rall(
    bpds::Vector{Float64}, 
    days::Vector{Float64}, 
    r_meters::Vector{Float64}, 
    S::Float64, T::Float64, rho::Float64;  # ✅ Separate parameters
    evaluation_days_from_start=nothing
)
```

### 2. **Coordinate Grid Confusion**
The monthly notebook had confusing coordinate handling with "flipped" grids:

**Before (Monthly - CONFUSING):**
```julia
# Call with lon_grid, lat_grid (swapped!)
pfield = pfieldcalc_all_rates(
    lon_grid, lat_grid, STRho, days, rates,  # ❌ Swapped coordinates
    inj[:well_lon], inj[:well_lat], "latlon",
    eval_days
)

# Inside function: "xGrid_km actually contains latitude values"
R_km = haversine_distance.(xGrid_km, yGrid_km, ywell_km, xwell_km)  # ❌ Confusing
```

**After (Monthly - CLEAR):**
```julia
# Call with lat_grid, lon_grid (natural order)
pfield = pfieldcalc_all_rates(
    lat_grid, lon_grid, S, T, rho, days, rates,  # ✅ Natural order
    inj[:well_lon], inj[:well_lat]; 
    evaluation_days_from_start=eval_days
)

# Inside function: straightforward distance calculation
dist_km = haversine_distance.(lat_grid, lon_grid, well_lat, well_lon)  # ✅ Clear
```

### 3. **Superposition Implementation Differences**
While both implementations were attempting to do superposition correctly, the monthly version had unnecessary complexity with pre-built arrays that made it harder to verify correctness.

## Changes Made

### Cell 3: Helper Functions
- ✅ Changed `pressureScenario_Rall` to use separate `S, T, rho` parameters instead of tuple
- ✅ Simplified superposition logic to match the working annual notebook
- ✅ Changed `pfieldcalc_all_rates` to use natural `lat_grid, lon_grid` parameter order
- ✅ Removed confusing coordinate "flipping" logic
- ✅ Used straightforward haversine distance calculation

### Cell 6: Run Julia Hydrology Model
- ✅ Removed `STRho` tuple creation
- ✅ Updated `run_year` function to accept separate `S, T, rho` parameters
- ✅ Fixed `pfieldcalc_all_rates` call to use natural coordinate order
- ✅ Changed to keyword argument syntax for `evaluation_days_from_start`

## Expected Results

After these fixes, the monthly regression notebook should now:

1. **Match the annual notebook's approach** - using the same proven logic
2. **Show minimal differences** near the well (RMSE < 0.1 psi, Max error < 1 psi)
3. **Have clearer, more maintainable code** - no coordinate confusion
4. **Properly implement Theis superposition** - matching MATLAB's monthly gold files

## Testing

To verify the fix works:

1. Open `fsp_hydrology_monthly_regression.ipynb` in Jupyter
2. Run all cells
3. Check the metrics table (Cell 7) - should show low RMSE values
4. Check the heatmaps (Cell 8) - difference plots should be mostly white/near-zero

## Why Differences Appeared Near the Well

The near-well differences were particularly noticeable because:

1. **Pressure gradients are steepest near the injection point** - small coordinate errors get magnified
2. **Coordinate system confusion** - if lat/lon were swapped, distance calculations would be wrong
3. **High sensitivity** - the Theis solution is very sensitive to distance near the source

The fixes ensure that:
- Distance calculations use the correct coordinate order
- The superposition logic matches the proven annual implementation
- Function signatures are consistent and clear

## References

- Working reference: `regression_fsp_hydrology_modeling.ipynb` (annual constant injection)
- Fixed notebook: `fsp_hydrology_monthly_regression.ipynb` (monthly FSP format)
- Gold data: `regression_tests/gold/MonthlyInjectionReference*.csv`




