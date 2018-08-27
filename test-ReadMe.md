#0.prepare
##----when reboot huawei machine for C test:----
>   	$sysctl -w vm.nr_hugepages=8192

##----when reboot huawei machine for JAVA test:----
>   	$nohup ./identify_server &

##----*.sh config----
>      add specified files in right dir

#1.for C test

 ##----compile----
>       dir: /dpdk_app 
>       use: 
>       	 $source *.sh (once) 
>       	 $make clean all(after source operation)
>   	to compile the project to get *.so
>   	
>   	then use:
>   		 $cp .so /hust/runtest/JNILib
>   	to copy *.so into the library

##----file need copied dir----
>   	/include : regs_info.h
>   	/example3/test : fpga_ddr_rw.c ; queue_ctrl.h ; ringbuffer.h

##----fpga program----
>   		$FpgaCmdEntry LF -S 0 -I [FPGA-IMAGE-ID]
>   		use:
>   			$fis fpga-image-id 
>   		to show the AEI image in OBS

##----software run----
>   	dir: /runtest/run.sh
>   	use:
>   		$./run.sh >./log/log.txt &
>   	to run the shell background and save the output message into the .txt file
>   	
>   	then use:
>   		$grep "map taskID" -c [file name]
>   		$grep "reduce taskID" -c [file name]
>   	to dubug whether the seeds are blocked in the RAM
>
>   	also can change following content in the run.sh:
>   		-Xmx30g ----> -Xmx15g
>   	if it runs faster to be blocked,that means seeds are blocked.
>   	
>   	use:
>   		$jps (to show the java job id which is running background)
>   		$kill -9 pid [pid] (to kill the background job)
>   	to kill the job
>   	
>   	note:
>   		when you are running 4_channel,change following content in run.sh:
>   			-nct x ----> -nct 2
>   		when you are running 12_channel,change following content in run.sh:
>   			-nct x ----> -nct 6
>   		when you are running 18_channel,change following content in run.sh:
>   			-nct x ----> -nct 9

#2.for JAVA test

##----JAVA-LOY(prepare pkg test)----
>   	in /example3/test : org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.c
>   		change JNI function :
>   			JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_pushData
>   	        (JNIEnv * env, jobject obj, jbyteArray write_buff){}
>   		make this function do nothing
>
>   	in /example3/test : org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.c
>   		change JNI function :
>               JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_resiveFromFPGA
>               (JNIEnv *env, jobject obj){
>   		make this function do nothing

##----api call----
>      in /example3/test : org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.c
>   		change JNI function :
>   			JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_pushData
>        (JNIEnv * env, jobject obj, jbyteArray write_buff){}
>   	    make this function always do 
>              writePkg(a); 
>           operation.
>
>      in /example3/test : queue_ctrl.h
>           change function :
>              void callback(unsigned int thread_id, unsigned int slot_id, rw_ddr_data rw_data, int rw_flag){}
>           delete 202 Row : updateDsTail();
>           make the DsTail not update

##----update regs_Info----
>      in /example3/test : org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.c
>   		change JNI function :
>   			JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_pushData
>        (JNIEnv * env, jobject obj, jbyteArray write_buff){}
>   		make this function do 
>              if (isWriteable){
>                  writePkg(a);
>              }
>              operation.
>
>      in /example3/test : queue_ctrl.h
>           change function :
>               void callback(unsigned int thread_id, unsigned int slot_id, rw_ddr_data rw_data, int rw_flag){}
>           add 202 Row : updateDsTail();
>           make the DsTail update

##----JAVA-LOY-HUAWEI(pkg test and optimization)----
>   	  do all