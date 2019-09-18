#!/usr/bin/env bash

# porting_test.sh

# check the system dependencies for pyprecag are installed.
# build a temporary virtual environment for python 2 and for
# python 3 and attempt installation in each.

if [[ "$1" == "2" ]] ; then
    TARGET_PYTHONS=("python2.7")
elif [[ "$1" == "3" ]] ; then
    TARGET_PYTHONS=("python3.6")
else
    TARGET_PYTHONS=("python2.7" "python3.6")
fi

set -x

function check_requirements {
    set -e
    gdal-config --version
    for P in "${TARGET_PYTHONS[@]}" ; do
        "${P}" --version
    done
    python3 -c 'import virtualenv'
    echo testing for python3.6-dev
    python3.6-config --configdir
    echo OK
    echo checking we are in a directory with setup.py
    test -r ./setup.py
    echo OK
    set +e

}

function install_gdal {
    _pip=$1
    "${_pip}" install GDAL=="$(gdal-config --version)" --global-option=build_ext --global-option=-I/usr/include/gdal

}

function install_pvtools {
    _pip=$1
    "${_pip}" install -e .[dev,test]
}


function run_tests {
    _python=$1
    timeout 150s time "${_python}" -m unittest discover -s pyprecag/tests -p "*test*.py" -v
}

function run_single_tests {
    # we have a test (or possibly setup) that is hanging. This is to try and identify which one
    _python=$1
    for _test in $(find pyprecag/tests/ -name 'test_*.py') ; do
        echo testing "${_test}"
        "${_python}" -m unittest "${_test}"
    done

}

check_requirements

TEMPDIR=$(mktemp -d)

for PYTHON in "${TARGET_PYTHONS[@]}" ; do
  python3 -m virtualenv -p "${PYTHON}" "${TEMPDIR}/${PYTHON}"
  _pip="${TEMPDIR}/${PYTHON}/bin/pip"
  ls -l "${_pip}"
  "${_pip}" --version
  install_gdal "${_pip}"
  install_pvtools "${_pip}"
  _python="${TEMPDIR}/${PYTHON}/bin/python"
  #run_tests "${_python}"
  run_single_tests "${_python}"
done

#find "${TEMPDIR}" -maxdepth 2
