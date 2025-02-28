cd /gpfs/wolf2/cades/cli185/scratch/ywo/E3SM/inputdata/lnd/clm2/surfdata_map

reorder_longxy() {
  local input_file=$1
  local temp_flipped="temp_flipped.nc"
  local temp_with_coord="temp_with_coord.nc"
  local temp_reordered="temp_reordered.nc"
  local output_file=$2

  # Convert the LONGXY from -180 to 180, to 0-360
  ncap2 -s 'where(LONGXY < 0) LONGXY=LONGXY+360' ${input_file} ${temp_flipped}

  # Create the lsmlon dimension
  ncap2 -O -s 'lsmlon[$lsmlon]=array(0,1,$lsmlon)' ${temp_flipped} ${temp_with_coord}

  # Reorder along the lsmlon dimension
  ncks -O --msa -d lsmlon,360,719 -d lsmlon,0,359 ${temp_with_coord} ${temp_reordered}

  # Delete the lsmlon dimension
  ncks -C -O -x -v lsmlon ${temp_reordered} ${output_file}

  # Delete temporary files
  rm ${temp_flipped} ${temp_with_coord} ${temp_reordered}
}

# concatenate the historical and SSP585 file
# (note the last and first years are both 2015)
ncks -d time,1, landuse.timeseries_0.5x0.5_ssp5_rcp85_simyr2015-2100_c240408.nc trimmed_second.nc
ncrcat landuse.timeseries_0.5x0.5_hist_simyr1850-2015_c240308.nc \
       trimmed_second.nc \
       -o landuse.timeseries_0.5x0.5_combined_simyr1850-2100_c240508.nc
rm trimmed_second.nc

reorder_longxy landuse.timeseries_0.5x0.5_combined_simyr1850-2100_c240508.nc landuse.timeseries_0.5x0.5_combined_simyr1850-2100_c240508_newlon.nc
reorder_longxy surfdata_0.5x0.5_simyr1850_c240308.nc surfdata_0.5x0.5_simyr1850_c240308_newlon.nc