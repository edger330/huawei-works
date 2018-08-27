HOMEDIR=`pwd`/..
myPath=${HOMEDIR}/JNILib
if [ -d "$myPath" ]; then
rm -rf $myPath
fi
mkdir ${HOMEDIR}/JNILib
#cd ${HOMEDIR}/JNILib
cd ../native
cd BatchSmithWaterman
make clean
make all
cp libBatchSmithWaterman.so ../../JNILib
cd ../SmithWaterman
make clean
make all
cp libSmithWaterman.so ../../JNILib
cd ../CPPPairHMM
make clean
make all
cp libCPPPairHMM.so ../../JNILib
cd ../FPGAPairHMM
make clean
make all
cp libFPGAPairHMM.so ../../JNILib
cd ../BaseRecalibrationCmethodFastMuTect2
make clean
make build
cp libBaseRecalibrationFastMuTect2.so ../../JNILib
cd ../BaseRecalibratorCmethod
make clean
make build
cp liborg_broadinstitute_gatk_utils_baq_BAQ.so ../../JNILib
cd ../BaseRecalibrationCmethodFastHaplotypeCaller
make clean
make build
cp libBaseRecalibrationFastHaplotypeCaller.so ../../JNILib
