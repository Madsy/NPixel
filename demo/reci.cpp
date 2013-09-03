#include <linealg.h>

int reci11(int val)
{
	int ret;
	asm volatile (
	     "xorl %%edx, %%edx;\n"
	     "movl $0x400000, %%eax;\n"
	     "idiv %0;\n"
	     "movl  %%eax, %1\n"
	     :"=r"(ret)
	     :"r"(val)
	     :"%eax", "%edx"
	);
	return ret;
}


int reci15(int val)
{
	int ret;
	asm volatile (
	     "xorl %%edx, %%edx;\n"
	     "movl $0x40000000, %%eax;\n"
	     "idiv %0;\n"
	     "movl  %%eax, %1\n"
	     :"=r"(ret)
	     :"r"(val)
	     :"%eax", "%edx"
	);
	return ret;
}

int reci8(int val)
{
	int ret;
	asm volatile (
	     "xorl %%edx, %%edx;\n"
	     "movl $0x100, %%eax;\n"
	     "idiv %0;\n"
	     "movl  %%eax, %1\n"
	     :"=r"(ret)
	     :"r"(val)
	     :"%eax", "%edx"
	);
	return ret;
}
