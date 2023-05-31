#run this in the directory you would like to save the data to

import urllib.request
import json


# Construct the URL. We'll use the jsonfilelist service
url = "https://archive.gemini.edu/jsonfilelist/"

# List the files for GN-2014A-DD-3 taken with GMOS-N from 2014-06-01 to 2014-07-01
url += "canonical/GN-2014A-DD-3/GMOS-N/20140601-20140701"

# Open the URL and fetch the JSON document text into a string
u = urllib.request.urlopen(url)
jsondoc = u.read()
u.close()

# Decode the JSON
files = json.loads(jsondoc)

# This is a list of dictionaries each containing info about a file
for f in files:
    print("Filename: %s" % f['filename'])
    print("-- file size: %d, data size: %d" % (f['file_size'], f['data_size']))

#downalod the data
for id_file, filex in enumerate(files):
    a_filename = filex['name']
    response = urllib.request.urlopen("http://archive.gemini.edu/file/" + a_filename)

    fits_bytes = response.read()
    response.close()

    with open(a_filename, 'wb') as f:
        f.write(fits_bytes)
    print("Downloaded file %d of %d" % (id_file + 1, len(files)))

#the reduced data is provided in the tables/ folder