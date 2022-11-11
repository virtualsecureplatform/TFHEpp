	.file	"spqlios-ifft-avx.s"
	.text
	.p2align 4
#if !__APPLE__
	.globl	ifft
	.type	ifft, @function
ifft:
#else
	.globl	_ifft
_ifft:
#endif

//typedef struct  {
//    uint64_t n;
//    double* trig_tables;
//} IFFT_PRECOMP;

/* void _ifft(const void *tables, double *real) */
	/* Save registers */
	pushq       %r10
	pushq       %r11
	pushq       %r12
	pushq       %r13
	pushq       %r14
	pushq       %rbx
	
	/* Permute registers for better variable names */
	movq        %rdi, %rax
	movq        %rsi, %rdi      /* rsi: base of the real data */
	
        //IFFT_PRECOMP* fft_tables = (IFFT_PRECOMP*) tables;
        //const int32_t n = fft_tables->n;
        //const double* trig_tables = fft_tables->trig_tables;

	/* Load struct FftTables fields */
	movq         0(%rax), %rdx  /* rdx: Size of FFT (a power of 2, must be at least 4) */
	movq         8(%rax), %r8   /* r8: Base address of trigonometric tables array */
	
        //int32_t ns4 = n/4;
        //double* are = c;    //size n/4 (x8 because doubles)
        //double* aim = c+ns4; //size n/4
	movq	%rdx, %r10
	shl	$1, %r10
	add	%r10, %rsi          /* rsi: base of the imaginary data */

	//    //multiply by omega^j
	//    for (int32_t j=0; j<ns4; j+=4) {
	//	const double* r0 = trig_tables+2*j;
	//	const double* r1 = r0+4;
	//	//(re*cos-im*sin) + i (im*cos+re*sin)
	//	double* d0 = are+j;
	//	double* d1 = aim+j;
	//	dotp4(tmp0,d0,r0); //re*cos
	//	dotp4(tmp1,d1,r0); //im*cos
	//	dotp4(tmp2,d0,r1); //re*sin
	//	dotp4(tmp3,d1,r1); //im*sin
	//	sub4(d0,tmp0,tmp3);
	//	add4(d1,tmp1,tmp2);
	//    }

	shr	$3, %r10  /* now, r10 is n/4 (the last iteration) */
	movq    $0, %rcx  /* Loop counter: Range [0, r10), step size 4 */
	movq	%r8, %r11 /* r11 is the trig table pointer, step size 64 */
firstloop:
        vmovapd (%rdi,%rcx,8), %zmm0 /* real */
        vmovapd (%rsi,%rcx,8), %zmm1 /* imag */
	vmovapd 0(%r11), %zmm2 /* cos */
	vmovapd 64(%r11), %zmm3 /* sin */
        vmulpd  %zmm0, %zmm2, %zmm4 /* re*cos */
        vmulpd  %zmm0, %zmm3, %zmm5 /* re*sin */
	vfnmadd231pd  %zmm1, %zmm3, %zmm4 /* re*cos - im*sin */
	vfmadd231pd  %zmm1, %zmm2, %zmm5 /* re*sin + im*cos */
	vmovapd %zmm4, (%rdi,%rcx,8)
	vmovapd %zmm5, (%rsi,%rcx,8)
        //next iteration
	leaq  128(%r11), %r11
	addq	$8,%rcx
	cmpq	%r10,%rcx
	jb	firstloop


/*	
    const double* cur_tt = trig_tables;
    for (int32_t nn=ns4; nn>=8; nn/=2) {
	int32_t halfnn = nn/2;
	cur_tt += 2*nn;
	for (int32_t block=0; block<ns4; block+=nn) {
	    for (int32_t off=0; off<halfnn; off+=4) {
		double* d00 = are + block + off;
		double* d01 = aim + block + off;
		double* d10 = are + block + halfnn + off;
		double* d11 = aim + block + halfnn + off;
		add4(tmp0,d00,d10); // re + re
		add4(tmp1,d01,d11); // im + im
		sub4(tmp2,d00,d10); // re - re
		sub4(tmp3,d01,d11); // im - im
		copy4(d00,tmp0);
		copy4(d01,tmp1);
		const double* r0 = cur_tt+2*off;
		const double* r1 = r0+4;
		dotp4(tmp0,tmp2,r0); //re*cos
		dotp4(tmp1,tmp3,r1); //im*sin
		sub4(d10,tmp0,tmp1);
		dotp4(tmp0,tmp2,r1); //re*sin
		dotp4(tmp1,tmp3,r0); //im*cos
		add4(d11,tmp0,tmp1);
	    }
	}
    }
*/
        /* r10 is still = n/4  (constant) */
        /* r8 is cur_tt (initially, base of trig table) */
	movq %r10,%r12 /* r12 (nn): outer loop counter from n/4 to 8 */
nnloop:
	movq %r12,%r13 
	shr  $1,%r13   /* r13 = halfnn */
	leaq (%r8,%r12,8),%r8  
	leaq (%r8,%r12,8),%r8 /* update cur_tt += nn*16 */ 
	movq $0,%r11   /* r11 (block) */
blockloop:
	leaq (%rdi,%r11,8),%rax /* are + block */
	leaq (%rsi,%r11,8),%rbx /* aim + block */
	leaq (%rax,%r13,8),%rcx /* are + block + halfnn */
	leaq (%rbx,%r13,8),%rdx /* aim + block + halfnn */
	movq $0,%r9   /* r9 (off) */
	movq %r8,%r14 /* r14 : cur_tt + 16*off */
offloop:
        vmovapd (%rax,%r9,8), %zmm0  /* re0 */
        vmovapd (%rbx,%r9,8), %zmm1  /* im0 */
        vmovapd (%rcx,%r9,8), %zmm2  /* re1 */
        vmovapd (%rdx,%r9,8), %zmm3  /* im1 */
	vaddpd %zmm0,%zmm2,%zmm4 /* re0+re1 */
	vaddpd %zmm1,%zmm3,%zmm5 /* im0+im1 */
	vsubpd %zmm2,%zmm0,%zmm6 /* re2=re0-re1 */
	vsubpd %zmm3,%zmm1,%zmm7 /* im2=im0-im1 */
	vmovapd %zmm4,(%rax,%r9,8)
        vmovapd %zmm5,(%rbx,%r9,8)
        vmovapd (%r14),%zmm8 /* cos */
	vmovapd 64(%r14),%zmm9 /* sin */
	vmulpd %zmm6,%zmm8,%zmm4 /* re2.cos */
	vfnmadd231pd %zmm7,%zmm9,%zmm4 /* re2.cos - im2.sin */
	vmulpd %zmm6,%zmm9,%zmm5 /* re2.sin */
	vfmadd231pd %zmm7,%zmm8,%zmm5 /* re2.sin + im2.cos */
	vmovapd %zmm4,(%rcx,%r9,8)
        vmovapd %zmm5,(%rdx,%r9,8)
        /* end of off loop */
       	leaq 128(%r14),%r14
        addq $8,%r9
	cmpq %r13,%r9
        jb offloop
	/* end of block loop */
	addq %r12,%r11
	cmpq %r10,%r11
	jb blockloop
	/* end of nn loop */
	movq %r13,%r12
	cmpq $16,%r12
	jae nnloop

// last iteration
	movq %r12,%r13 
	shr  $1,%r13   /* r13 = halfnn */
	leaq (%r8,%r12,8),%r8  
	leaq (%r8,%r12,8),%r8 /* update cur_tt += nn*16 */ 
	movq $0,%r11   /* r11 (block) */
blockloop2:
	leaq (%rdi,%r11,8),%rax /* are + block */
	leaq (%rsi,%r11,8),%rbx /* aim + block */
	leaq (%rax,%r13,8),%rcx /* are + block + halfnn */
	leaq (%rbx,%r13,8),%rdx /* aim + block + halfnn */
	movq %r8,%r14 /* r14 : cur_tt + 16*off */
        vmovapd (%rax), %ymm0  /* re0 */
        vmovapd (%rbx), %ymm1  /* im0 */
        vmovapd (%rcx), %ymm2  /* re1 */
        vmovapd (%rdx), %ymm3  /* im1 */
	vaddpd %ymm0,%ymm2,%ymm4 /* re0+re1 */
	vaddpd %ymm1,%ymm3,%ymm5 /* im0+im1 */
	vsubpd %ymm2,%ymm0,%ymm6 /* re2=re0-re1 */
	vsubpd %ymm3,%ymm1,%ymm7 /* im2=im0-im1 */
	vmovapd %ymm4,(%rax)
        vmovapd %ymm5,(%rbx)
        vmovapd (%r14),%ymm8 /* cos */
	vmovapd 32(%r14),%ymm9 /* sin */
	vmulpd %ymm6,%ymm8,%ymm4 /* re2.cos */
	vfnmadd231pd %ymm7,%ymm9,%ymm4 /* re2.cos - im2.sin */
	vmulpd %ymm6,%ymm9,%ymm5 /* re2.sin */
	vfmadd231pd %ymm7,%ymm8,%ymm5 /* re2.sin + im2.cos */
	vmovapd %ymm4,(%rcx)
        vmovapd %ymm5,(%rdx)

	/* end of block loop */
	addq %r12,%r11
	cmpq %r10,%r11
	jb blockloop2


/*
    //size 4 loop
    {
	for (int32_t block=0; block<ns4; block+=4) {
	    double* d0 = are+block;
	    double* d1 = aim+block;
	    tmp0[0]=d0[0];
	    tmp0[1]=d0[1];
	    tmp0[2]=d0[0];
	    tmp0[3]=-d1[1];
	    tmp1[0]=d0[2];
	    tmp1[1]=d0[3];
	    tmp1[2]=-d0[2];
	    tmp1[3]=d1[3];
	    tmp2[0]=d1[0];
	    tmp2[1]=d1[1];
	    tmp2[2]=d1[0];
	    tmp2[3]=d0[1];
	    tmp3[0]=d1[2];
	    tmp3[1]=d1[3];
	    tmp3[2]=-d1[2];
	    tmp3[3]=-d0[3];
	    add4(d0,tmp0,tmp1);
	    add4(d1,tmp2,tmp3);
	}
    }
*/
    	/* r10 is still = n/4  (constant) */
	vmovapd     size4negation0(%rip), %zmm15
	vmovapd     size4negation1(%rip), %zmm14
	vmovapd     size4negation2(%rip), %zmm13
	vmovapd     size4negation3(%rip), %zmm12
	vmovapd     permutex1(%rip), %zmm11
	vmovapd     permutex2(%rip), %zmm10
	movq $0,%rax /* rax (block) */
	movq %rdi,%r11 /* r11 (are+block) */
	movq %rsi,%r12 /* r12 (aim+block) */
size4loop:
	vmovapd (%r11),%zmm0 /* r0 r1 r2 r3 */
	vmovapd (%r12),%zmm1 /* i0 i1 i2 i3 */

	# vshufpd $10,%zmm1,%zmm0,%zmm2 /* r0 i1 r2 i3 */
	# vshufpd $10,%zmm0,%zmm1,%zmm3 /* i0 r1 i2 r3 */
	# vperm2f128 $32,%zmm2,%zmm0,%zmm4 /* r0 r1 r0 i1 */
	# vperm2f128 $49,%zmm2,%zmm0,%zmm5 /* r2 r3 r2 i3 */
	# vperm2f128 $32,%zmm3,%zmm1,%zmm6 /* i0 i1 i0 r1 */
	# vperm2f128 $49,%zmm3,%zmm1,%zmm7 /* i2 i3 i2 r3 */
	vmovapd %zmm11, %zmm4
	vmovapd %zmm11, %zmm6
	vmovapd %zmm10, %zmm5
	vmovapd %zmm10, %zmm7
	vpermi2pd %zmm1, %zmm0, %zmm4
	vpermi2pd %zmm0, %zmm1, %zmm6
	vpermi2pd %zmm1, %zmm0, %zmm5
	vpermi2pd %zmm0, %zmm1, %zmm7

	vmulpd	%zmm4,%zmm15,%zmm4 /* r0 r1 r0 -i1 */
	vfmadd231pd	%zmm5,%zmm14,%zmm4 /* (r0 r1 r0 -i1) + (r2 r3 -r2 i3) */
	vfmadd231pd	%zmm7,%zmm13,%zmm6 /* (i0 i1 i0 r1) + (i2 i3 -i2 -r3) */
	vmovapd %zmm4,(%r11) 
	vmovapd %zmm6,(%r12)
        /* end of loop */
        leaq 64(%r11),%r11
        leaq 64(%r12),%r12
        addq $8,%rax
	cmpq %r10,%rax
	jb size4loop


/*
    //size 2
    {
	for (int32_t block=0; block<ns4; block+=4) {
	    double* d0 = are+block;
	    double* d1 = aim+block;
	    tmp0[0]=d0[0];
	    tmp0[1]=d0[0];
	    tmp0[2]=d0[2];
	    tmp0[3]=d0[2];
	    tmp1[0]=d0[1];
	    tmp1[1]=-d0[1];
	    tmp1[2]=d0[3];
	    tmp1[3]=-d0[3];
	    add4(d0,tmp0,tmp1);
	    tmp0[0]=d1[0];
	    tmp0[1]=d1[0];
	    tmp0[2]=d1[2];
	    tmp0[3]=d1[2];
	    tmp1[0]=d1[1];
	    tmp1[1]=-d1[1];
	    tmp1[2]=d1[3];
	    tmp1[3]=-d1[3];
	    add4(d1,tmp0,tmp1);
	}
    }
}
*/
	movq $0,%rax /* rax (block) */
	movq %rdi,%r11 /* r11 (are+block) */
	movq %rsi,%r12 /* r12 (aim+block) */
size2loop:
	vmovapd (%r11),%zmm0 /* r0 r1 r2 r3 */
	vmovapd (%r12),%zmm1 /* i0 i1 i2 i3 */
	vshufpd $0,%zmm0,%zmm0,%zmm2 /* r0 r0 r2 r2 */
	vshufpd $255,%zmm0,%zmm0,%zmm3 /* r1 r1 r3 r3 */
	vshufpd $0,%zmm1,%zmm1,%zmm4 /* i0 i0 i2 i2 */
	vshufpd $255,%zmm1,%zmm1,%zmm5 /* i1 i1 i3 i3 */
	vfmadd231pd %zmm3,%zmm12,%zmm2 /* (r0 r0 r2 r2) + (r1 -r1 r3 -r3) */
	vfmadd231pd %zmm5,%zmm12,%zmm4 /* (i0 i0 i2 i2) + (i1 -i1 i3 -i3) */ 
	vmovapd %zmm2,(%r11)
	vmovapd %zmm4,(%r12)
	/* end of loop */
        leaq 64(%r11),%r11
        leaq 64(%r12),%r12
        addq $8,%rax
	cmpq %r10,%rax
	jb size2loop


	/* Restore registers */
end:
	vzeroall
	popq        %rbx
	popq        %r14
	popq        %r13
	popq        %r12
	popq        %r11
	popq        %r10
	retq


/* Constants for YMM */
.balign 64
size4negation0: .double +1.0, +1.0, +1.0, -1.0, +1.0, +1.0, +1.0, -1.0
size4negation1: .double +1.0, +1.0, -1.0, +1.0, +1.0, +1.0, -1.0, +1.0
size4negation2: .double +1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0
size4negation3: .double +1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +1.0, -1.0
permutex1: .quad 0, 1, 0, 8+1, 4, 5, 4, 8+5 /* zmm11 */
permutex2: .quad 2, 3, 2, 8+3, 6, 7, 6, 8+7 /* zmm10 */


#if !__APPLE__
	.size	ifft, .-ifft
#endif

