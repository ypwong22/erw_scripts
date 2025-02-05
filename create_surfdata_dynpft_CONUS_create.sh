create=0
submit=0
clean=1

if [[ $create == 1 ]]; then
  for chunk in {1..23}; do # SHOULD count until 190 instead of 4
    sed "s/REPLACE/${chunk}/g" create_surfdata_dynpft_CONUS.py > temp_create_surfdata_dynpft_CONUS_${chunk}.py
    sed "s/REPLACE/${chunk}/g" create_surfdata_dynpft_CONUS.sh > temp_create_surfdata_dynpft_CONUS_${chunk}.sh
  done
fi

if [[ $submit == 1 ]]; then
  for chunk in {1..23}; do
    sbatch temp_create_surfdata_dynpft_CONUS_${chunk}.sh
  done
fi

if [[ $clean == 1 ]]; then
  for chunk in {1..23}; do
    rm temp_create_surfdata_dynpft_CONUS_${chunk}.py
    rm temp_create_surfdata_dynpft_CONUS_${chunk}.sh
  done
fi
