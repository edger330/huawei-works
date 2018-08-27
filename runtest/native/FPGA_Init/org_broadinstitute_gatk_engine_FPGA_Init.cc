#include "org_broadinstitute_gatk_engine_FPGA_Init.h"
#include <iostream>
#include <sys/socket.h>
#include <memory.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <zconf.h>
#include <cstdio>
#include <sys/types.h>
#include <cstdlib>
#include <cstring>
#include <netinet/in.h>
#include <cerrno>

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

JNIEXPORT void JNICALL Java_org_broadinstitute_gatk_engine_FPGA_1Init_init(JNIEnv *, jobject) {
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

    if (!correct) exit(1);
}