#!/bin/bash

set -e

# Created 2024-09-23 11:33:02

CASEDIR="/glade/u/home/geraint/NilsScheme/NilsFSAdapt"

/glade/work/bitz/cesm2_3_beta17wave_06242024/cime/scripts/create_newcase --case "${CASEDIR}" --compset 2000_DATM%JRA-1p4-2018_SLND_CICE_DOCN%SOM_DROF%JRA-1p4-2018_SGLC_WW3DEV_SESP --res TL319_t232_wg37 --run-unsupported --project UWAS0136

cd "${CASEDIR}"

./case.setup

./pelayout

./case.build

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build --clean ice

./case.build

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.build

./case.build

./case.build --clean ice

./case.build

./case.build

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.build

./case.build

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean-all

./case.setup --clean

./case.setup

./pelayout

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.build

./case.submit

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build --clean ice

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build --clean ice

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build

./case.build --clean ice

./case.build

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

./case.build --clean ice

./case.build

./case.submit

