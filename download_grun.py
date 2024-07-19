"""
Ghiggi, G., Humphrey, V., Seneviratne, S. I., & Gudmundsson, L. (2021). G-RUN ENSEMBLE: A multi-forcing observation-based global runoff reanalysis. Water Resources Research, 57(5), e2020WR028787. https://doi.org/10.1029/2020WR028787

- G-RUN_ENSEMBLE_MMM.nc covers the time period from 1902 to 2019 and provide the median of the G-RUN ENSEMBLE members. If you want to rely on one single estimate this is likely the file you are interested in.

- G-RUN_ENSEMBLE_MEMBERS.zip contains ensemble mean reconstructions for 21 different atmospheric forcing datasets. The time range depends on the considered forcing.

- Each remaining file called G-RUN_ENSEMBLE_*.zip (where * denotes the acronym of the atmospheric forcing dataset used to force the model), contains 25 runoff reconstructions obtained by training models on different subsets of the available runoff observations.
"""
import requests
import os
from tqdm import tqdm
import xarray as xr


# Converting the file content from JSON string to a Python data structure
def fetch_json(url):
    try:
        response = requests.get(url)
        response.raise_for_status()  # This will raise an exception for HTTP errors
        return response.json()  # Parse and return the JSON data
    except requests.RequestException as e:
        print(f"An error occurred: {e}")
        return None


def download_files(url):
    json_data = fetch_json(url)

    for file in json_data['files']:
        download_url = file['download_url']
        filename = os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'GRUN', file['name'])

        if os.path.exists(filename):
            continue

        response = requests.get(download_url, stream=True)
        if response.status_code == 200:
            with open(filename, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"File downloaded: {file['name']}")
        else:
            print(f"Failed to download file, status code: {response.status_code}")


if __name__ == '__main__':
    url_grun = "https://api.figshare.com/v2/articles/12794075"
    download_files(url_grun)

    # subset to CONUS
    grun = xr.open_dataset(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'GRUN',
                                        'G-RUN_ENSEMBLE_MMM.nc'))

    time_slice = slice('1990-01-01', '2009-12-31')  # Replace with your desired time range
    bounding_box = {'lon': slice(-125.25, -66.5), 'lat': slice(23.25, 54.7)}
    runoff_obs = grun['Runoff'].sel(time=time_slice, X=bounding_box['lon'], 
                                    Y=bounding_box['lat']).load()
    runoff_obs = runoff_obs.rename({'X': 'lon', 'Y': 'lat'})
    runoff_obs.to_netcdf(os.path.join(os.environ['PROJDIR'], 'ERW_LDRD', 'data', 'GRUN',
                                        'G-RUN_ENSEMBLE_MMM_CONUS.nc'))
    grun.close()