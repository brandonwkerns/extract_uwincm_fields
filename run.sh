#!/bin/sh

python extract_uwincm_fields.py /home/orca/from_miami/uwin1/milan/output/alaska/awo_2014110700_fnl_4km_nudging
cp alaska_awo_2014110700_fnl_4km_nudging.nc4 /home/orca/bkerns/public_html/sharedata/alaska_output

exit 0

