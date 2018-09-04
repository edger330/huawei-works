#!/bin/bash

gatk=JavaFpgaTest_main.jar

lib_path=JNILib

time java \
-Djava.library.path=$lib_path -jar $gatk 
