#ifndef SRC_DUPLEX_QUEUE_DUPLEX_QUEUE_H_
#define SRC_DUPLEX_QUEUE_DUPLEX_QUEUE_H_

#define DS_PSIZE_EXP	17	// package size 128KiB
#define DS_PCAP_EXP		11	// package capacity 2048 (2047 actual available)

#define US_PSIZE_EXP	17	// package size 128KiB
#define US_PCAP_EXP		11	// package capacity 2048 (2047 actual available)

#define DS_PSIZE		(1L << DS_PSIZE_EXP)
#define DS_PCAP			((1L << DS_PCAP_EXP) - 1)
#define DS_ADDR_LIMIT	(DS_BASE_ADDR + (1L << (DS_PSIZE_EXP + DS_PCAP_EXP)))

#define US_PSIZE		(1L << US_PSIZE_EXP)
#define US_PCAP			((1L << US_PCAP_EXP) - 1)
#define US_ADDR_LIMIT	(US_BASE_ADDR + (1L << (US_PSIZE_EXP + US_PCAP_EXP)))
#define INTS_NUM	(1024 * 1024 / 4)

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>

#define DS_BASE_ADDR	0x00000000
#define US_BASE_ADDR	0x40000000

#define MAP_SIZE (32*1024UL)
#define MAP_MASK (MAP_SIZE - 1)

#define DS_TAIL_REG		0
//#define DS_HEAD_REG     1
#define DS_CNT_REG    	2
//#define US_TAIL_REG	  4
#define US_HEAD_REG		5
#define US_CNT_REG		6

void *regBase;
int regFd;
int dmaH2C, dmaC2H;

uint dsTail = DS_BASE_ADDR;
uint usHead = US_BASE_ADDR;

void writeReg(uint regNum, uint data)
{
    *(regNum + (uint *)regBase) = data;
}
uint readReg(uint regNum)
{
    return *(regNum + (uint *)regBase);
}
int OpenQueueCtrl()
{
    regFd = open("/dev/xdma0_user", O_RDWR | O_SYNC);
    regBase = mmap(0, MAP_SIZE, PROT_READ | PROT_WRITE, MAP_SHARED, regFd, 0);
    writeReg(DS_TAIL_REG, dsTail);
    writeReg(US_HEAD_REG, usHead);

    dmaH2C = open("/dev/xdma0_h2c_0", O_RDWR);
    dmaC2H = open("/dev/xdma0_c2h_0", O_RDWR | O_NONBLOCK);

    if(regFd >= 0 && dmaH2C >= 0 && dmaC2H >= 0)
    	return 1;
    else
    	return 0;
}
void CloseQueueCtrl()
{
	if(regFd >= 0)
		close(regFd);
	if(dmaH2C >= 0)
		close(dmaH2C);
	if(dmaC2H >= 0)
		close(dmaC2H);
}
int GetDsPkgCnt()
{
	return readReg(DS_CNT_REG);
}
int IsDsWriteable(int pkgCnt)
{
	int cnt = GetDsPkgCnt();
	return cnt + pkgCnt <= DS_PCAP;
}
void updateDsTail(int pkgCnt)
{
	for(; pkgCnt > 0; pkgCnt--)
	{
		dsTail += DS_PSIZE;
		if(dsTail >= DS_ADDR_LIMIT)
			dsTail = DS_BASE_ADDR;
	}
	writeReg(DS_TAIL_REG, dsTail);
}
void WriteDsPkg(void *buf, uint pkgCnt)
{
    off_t off = lseek(dmaH2C, dsTail, SEEK_SET);
    int rc = write(dmaH2C, buf, DS_PSIZE * pkgCnt);
    updateDsTail(pkgCnt);
}

int GetUsPkgCnt()
{
	return readReg(US_CNT_REG);
}
int IsUsReadable(int pkgCnt)
{
	int cnt = GetUsPkgCnt();
	return cnt >= pkgCnt;
}
void updateUsHead(int pkgCnt)
{
	for(; pkgCnt > 0; pkgCnt--)
	{
		usHead += US_PSIZE;
		if(usHead >= US_ADDR_LIMIT)
			usHead = US_BASE_ADDR;
	}
	writeReg(US_HEAD_REG, usHead);
}
void ReadUsPkg(void *buf, uint pkgCnt)
{
	off_t off = lseek(dmaC2H, usHead, SEEK_SET);
	int rc = read(dmaC2H, buf, US_PSIZE * pkgCnt);
	updateUsHead(pkgCnt);
}

#endif