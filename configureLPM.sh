#!/bin/bash

LPM_ROOT=/Users/pbosler/lpmCxx

rm -rf CMakeCache.txt CMakeFiles/

cmake \
-D CMAKE_INSTALL_PREFIX=/Users/pbosler/Desktop/lpmCxx/lpmInstall \
$LPM_ROOT

