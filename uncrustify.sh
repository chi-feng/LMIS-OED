#!/bin/bash

for f in src/*.cpp
do
 uncrustify -c uncrustify.cfg --replace $f
done

for f in src/*.h
do
 uncrustify -c uncrustify.cfg --replace $f
done

for f in src/pce/*.cpp
do
 uncrustify -c uncrustify.cfg --replace $f
done

for f in src/pce/*.h
do
 uncrustify -c uncrustify.cfg --replace $f
done

for f in src/merr/*.cpp
do
 uncrustify -c uncrustify.cfg --replace $f
done

for f in src/merr/*.h
do
 uncrustify -c uncrustify.cfg --replace $f
done

find src/*.unc-backup* | xargs rm -rvf
find src/pce/*.unc-backup* | xargs rm -rvf
find src/merr/*.unc-backup* | xargs rm -rvf