cmd_malloc_heap.o = gcc -Wp,-MD,./.malloc_heap.o.d.tmp -m64 -pthread -fPIC  -march=native -DRTE_MACHINE_CPUFLAG_SSE -DRTE_MACHINE_CPUFLAG_SSE2 -DRTE_MACHINE_CPUFLAG_SSE3 -DRTE_MACHINE_CPUFLAG_SSSE3 -DRTE_MACHINE_CPUFLAG_SSE4_1 -DRTE_MACHINE_CPUFLAG_SSE4_2 -DRTE_MACHINE_CPUFLAG_AES -DRTE_MACHINE_CPUFLAG_PCLMULQDQ -DRTE_MACHINE_CPUFLAG_AVX -DRTE_MACHINE_CPUFLAG_RDRAND -DRTE_MACHINE_CPUFLAG_FSGSBASE -DRTE_MACHINE_CPUFLAG_F16C -DRTE_MACHINE_CPUFLAG_AVX2  -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/build/include -include /home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/build/include/rte_config.h -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_eal/linuxapp/eal/include -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_eal/common -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_eal/common/include -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_ring -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_mempool -I/home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_ivshmem -W -Wall -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations -Wold-style-definition -Wpointer-arith -Wcast-align -Wnested-externs -Wcast-qual -Wformat-nonliteral -Wformat-security -Wundef -Wwrite-strings -O3   -o malloc_heap.o -c /home/hust/zzt/poc-fpga/fp1/software/userspace/dpdk_src/dpdk-16.04/lib/librte_eal/common/malloc_heap.c 
