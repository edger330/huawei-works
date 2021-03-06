/* DO NOT EDIT THIS FILE - it is machine generated */
#include <iostream>
#include <jni.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <stdint.h>
#include <sys/socket.h>
#include "utils.h"
#include "queue_ctrl.h"
#include <memory.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <zconf.h>



/* Header for class org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM */

#include "org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM.h"
/*
 * Class:     org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM
 * Method:    initFPGA
 * Signature: ()V
 */
#define _write_size 128*1024
#define _read_size 128*1024
#define S_PORT 20006
#define BUFFER_SIZE 1024
static int dev_fd=-1;
char write_buff[_write_size];
char read_buff[_read_size];
#define BUFFER_SIZE 16

using namespace std;

void get_randchar(unsigned char *rand_char) {
    srand((unsigned) (time(NULL)));
    for (int i = 0; i < BUFFER_SIZE; i++) {
        rand_char[i] = (unsigned char) (rand() % 256);
    }
}

void and_char(const unsigned char *rand_char, int *correctKey) {
    unsigned char pass[BUFFER_SIZE] = {
            127, 255, 203, 69,
            198, 123, 147, 209,
            75, 43, 252, 19,
            221, 87, 139, 55};

    for (int i = 0; i < BUFFER_SIZE; i++) {
        correctKey[(i / 4)] = (correctKey[(i / 4)] << 8) | (rand_char[i] & pass[i]);
    }
}

void init() {
    const char *server_ip = getenv("YT_CONNECT_IP");
    const char *server_port = getenv("YT_CONNECT_PORT");


    if (server_ip == NULL || server_port == NULL) {
        cout << "Please set the correct environment variable" << endl;
        exit(1);
    }

    int sock_cli = socket(AF_INET, SOCK_STREAM, 0);

    struct sockaddr_in servaddr;
    memset(&servaddr, 0, sizeof(servaddr));
    servaddr.sin_family = AF_INET;
    servaddr.sin_port = htons(atoi(server_port));
    servaddr.sin_addr.s_addr = inet_addr(server_ip);

    if (connect(sock_cli, (struct sockaddr *) &servaddr, sizeof(servaddr)) < 0) {
        perror("connect");
        exit(1);
    }

    unsigned char sendbuf[BUFFER_SIZE];
    get_randchar(sendbuf);

    int recvbuf[BUFFER_SIZE / 4];

    send(sock_cli, sendbuf, BUFFER_SIZE, 0); ///发送

    recv(sock_cli, recvbuf, sizeof(recvbuf), 0); ///接收

    int correctKey[BUFFER_SIZE / 4];

    and_char((unsigned char *) sendbuf, correctKey);

    bool correct = true;
    for (int i = 0; i < BUFFER_SIZE / 4; i++) {
        if (correctKey[i]!= recvbuf[i]) correct = false;
    }

    if (!correct){
        printf("error\n");
        exit(1);
    }
}

JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_initFPGA
  (JNIEnv *env, jobject obj){
     init();
     int open_statue = OpenQueueCtrl();
     if(!open_statue){
        printf("Open xdma device failed!\n");
        exit(-1);
     }
     printf("Open xdma deivce success!\n");

}

/*
 * Class:     org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM
 * Method:    doneFPGA
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_doneFPGA
        (JNIEnv *, jobject){
    CloseQueueCtrl();
}

int data_checkHex(uint8_t* buff){
    FILE * fp = fopen("./testDataFPGA.txt","a");
        uint8_t *dataPos = buff;
        dataPos +=15;
        uint8_t totalPairNum = XBUF_TO_UINT08(dataPos); dataPos += 1;
                fprintf(fp,"000000000000000000000000000000%02x",totalPairNum);
                for(uint8_t i = 0;i<totalPairNum;i++){
                    fprintf(fp,"\n");
                    dataPos+=14;
                    uint16_t PairLen = (XBUF_TO_UINT16(dataPos)*16); dataPos += 2;
                    int Len = (int)PairLen;
                    fprintf(fp,"0000000000000000000000000000%04x\n",PairLen);
                    for(int i=0;i<(Len/4);i++){
                        /*if((i%16)==0){
                            fprintf(fp,"\n");
                        }*/
                        uint32_t tempData = XBUF_TO_UINT32(dataPos);dataPos+=4;
                        fprintf(fp,"%08x ",tempData);
                    }
                }
                //fprintf(fp,"\nready to send data to FPGA\n");
                fprintf(fp,"\n");
                fprintf(fp,"\n");
               	fclose(fp);
                return 1;
}
int data_print(uint8_t* buff){
    FILE * fp = fopen("./testDataFPGA.txt","a");
    uint8_t *dataPos = buff;
    fprintf(fp,"package start!\n");
    for(int i=0;i<128*1024;i++){
        uint8_t temp = XBUF_TO_UINT08(dataPos); dataPos++;
        fprintf(fp,"%d\n",temp);
    }
    fprintf(fp,"\n");
}

int data_check(uint8_t* buff){
    FILE * fp = fopen("./testDataFPGA.txt","a");
    //    fprintf(fp,"===============================2=============================\n");
    //    fprintf(fp,"before write\n");


    /*
            uint32_t i = 0;
            uint8_t *dataPos = buff;
            uint32_t taskId = XBUF_TO_UINT32(dataPos); dataPos += 4;
            uint16_t sampleId = XBUF_TO_UINT16(dataPos); dataPos += 2;
            uint32_t subTaskId = XBUF_TO_UINT32(dataPos); dataPos += 4;
            dataPos += 14; // nc

           fprintf(fp ,"taskId:%u sampleId:%u subTaskId:%u\n", taskId, sampleId, subTaskId);*/
    //--------------------------------------------------------------------------------------//
            uint8_t *dataPos = buff;
            dataPos+=15;
            uint8_t totalPairNum = XBUF_TO_UINT08(dataPos); dataPos += 1;
            for(uint8_t i = 0;i<totalPairNum;i++){
                dataPos+=14;
                uint16_t PairLen = (XBUF_TO_UINT16(dataPos)*16); dataPos += 2;//
                uint16_t taskId = XBUF_TO_UINT16(dataPos); dataPos += 2;
                uint8_t sampleId = XBUF_TO_UINT08(dataPos); dataPos += 1;
                uint8_t hapIdx = XBUF_TO_UINT08(dataPos); dataPos += 1;
                uint16_t readIdx = XBUF_TO_UINT16(dataPos); dataPos += 12;//
                uint16_t hapLen = XBUF_TO_UINT16(dataPos); dataPos += 2;
                uint16_t rsLen = XBUF_TO_UINT16(dataPos); dataPos += 2;
                dataPos += (PairLen-22);
                fprintf(fp ,"PairLen:%u taskId:%u sampleId:%u hapIdx:%u readIdx:%u hapLen:%u readLen:%u\n",PairLen, taskId, sampleId, hapIdx,readIdx,hapLen,rsLen);
            }
	fclose(fp);

        return 1;
}


/*
 * Class:     org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM
 * Method:    sendDataToFPGA
 * Signature: ([B)V
 */
JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_sendDataToFPGA
  (JNIEnv * env, jobject obj, jbyteArray write_buff){
/*    jboolean is_copy = JNI_FALSE;
    jbyte * buffData = env->GetByteArrayElements(write_buff,&is_copy);
	char * a = (char*)buffData;
	uint8_t *tmp=(uint8_t *)a;
	data_check(tmp);
	wr
		fprintf(fp,"\n");ite(dev_fd,a,_write_size);
	env->ReleaseByteArrayElements(write_buff,buffData,JNI_COMMIT);
*/
//-----------------------------------------------------------------------//
	jbyte buffData[_write_size];
	env->GetByteArrayRegion(write_buff,0,_write_size,buffData);
	char *a = (char*)buffData;
	uint8_t *tmp = (uint8_t *)a;
	//data_print(tmp);
	//write(dev_fd,a,_write_size);
	while(true){
	    if(IsDsWriteable(1)){
	        WriteDsPkg(a,1);
	        break;
	    }
	}
//-----------------------------------------------------------------------//
}

/*
 * Class:     org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM
 * Method:    resiveFromFPGA
 * Signature: ()[B
 */
JNIEXPORT jbyteArray JNICALL Java_org_broadinstitute_gatk_utils_pairhmm_FPGAPairHMM_resiveFromFPGA
  (JNIEnv *env, jobject obj){
    memset(read_buff,0,_read_size);
    //read(dev_fd,read_buff,_read_size);
    while(true){
    	    if(IsUsReadable(1)){
    	        ReadUsPkg(read_buff,1);
    	        break;
    	    }
    	}
    jbyteArray result = env->NewByteArray(_read_size);
/*    FILE * fp = fopen("./testResultFPGA.txt","a");
//    fprintf(fp,"============================================================\n");
    fprintf(fp,"read from fpga\n");
   for(int i=0;i<16*19200;i++){
	if(i%16==0)
		fprintf(fp,"\n");
        fprintf(fp,"%02x ",(unsigned char)read_buff[i]);
    }
    fclose(fp);*/
    env->SetByteArrayRegion(result,0,_read_size,(jbyte*)read_buff);
    return result;
}
