	.file	"spqlios-fft-avx2.s"
#if !__APPLE__
	.section .note.GNU-stack,"",%progbits
#endif
	.text
	.p2align 4
#if !__APPLE__
	.globl	fft
	.type	fft, @function
fft:
#else
	.globl	_fft
_fft:
#endif
//c has size n/2
//void fft(const void* tables, double* c) {

//    FFT_PRECOMP* fft_tables = (FFT_PRECOMP*) tables;
//    const int32_t n = fft_tables->n;
//    const double* trig_tables = fft_tables->trig_tables;

	/* Save registers */
	pushq       %r10
	pushq       %r11
	pushq       %r12
	pushq       %r13
	pushq       %r14
	pushq       %rbx
	
	/* Permute registers for better variable names */
	movq        %rdi, %rax
	movq        %rsi, %rdi      /* rdi: base of the real data CONSTANT */
	
	/* Load struct FftTables fields */
	movq         0(%rax), %rdx  /* rdx: n (logical Size of _fft  = a power of 2, must be at least 4) */
	movq         8(%rax), %r8   /* r8: Base address of trigonometric tables array (CONSTANT) */
	
//    int32_t ns4 = n/4;
//    double* pre = c;     //size n/4
//    double* pim = c+ns4; //size n/4
	movq	%rdx, %r9
	shr	$2,%r9              /* r9: ns4 CONSTANT */
	leaq	(%rdi,%r9,8),%rsi   /* rsi: base of imaginary data CONSTANT */

//    //size 2
//    {
//	//[1  1]
//	//[1 -1]
//	//     [1  1]
//	//     [1 -1]
//	for (int32_t block=0; block<ns4; block+=4) {
//	    double* d0 = pre+block;
//	    double* d1 = pim+block;
//	    tmp0[0]=d0[0];
//	    tmp0[1]=d0[0];
//	    tmp0[2]=d0[2];
//	    tmp0[3]=d0[2];
//	    tmp1[0]=d0[1];
//	    tmp1[1]=-d0[1];
//	    tmp1[2]=d0[3];
//	    tmp1[3]=-d0[3];
//	    add4(d0,tmp0,tmp1);
//	    tmp0[0]=d1[0];
//	    tmp0[1]=d1[0];
//	    tmp0[2]=d1[2];
//	    tmp0[3]=d1[2];
//	    tmp1[0]=d1[1];
//	    tmp1[1]=-d1[1];
//	    tmp1[2]=d1[3];
//	    tmp1[3]=-d1[3];
//	    add4(d1,tmp0,tmp1);
//	}
//    }
	vmovapd     size4negation0(%rip), %zmm15
	vmovapd     size4negation1(%rip), %zmm14
	vmovapd     size4negation2(%rip), %zmm13
	vmovapd     size4negation3(%rip), %zmm12
	vmovapd     permutex1(%rip), %zmm11
	
	movq	$0,%rax	/* rax: block */
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize2loop:
	vmovapd (%r10),%zmm0 /* r0 r1 r2 r3 r4 r5 r6 r7 */
	vmovapd (%r11),%zmm1 /* i0 i1 i2 i3 i4 i5 i6 i7 */
	vshufpd $0,%zmm0,%zmm0,%zmm2  /* r0 r0 r2 r2 r4 r4 r6 r6 */
	vshufpd $255,%zmm0,%zmm0,%zmm3 /* r1 r1 r3 r3 r5 r5 r7 r7 */
	vshufpd $0,%zmm1,%zmm1,%zmm4  /* i0 i0 i2 i2 i4 i4 i6 i6 */
	vshufpd $255,%zmm1,%zmm1,%zmm5 /* i1 i1 i3 i3 i5 i5 i7 i7 */
	vfmadd231pd %zmm3,%zmm12,%zmm2 /* (r0 r0 r2 r2 r4 r4 r6 r6) + (r1 -r1 r3 -r3 r5 -r5 r7 -r7) */
	vfmadd231pd %zmm5,%zmm12,%zmm4 /* (i0 i0 i2 i2 i4 i4 i6 i6) + (i1 -i1 i3 -i3 i5 -i5 i7 -i7) */
	vmovapd %zmm2,(%r10)
	vmovapd %zmm4,(%r11)
	/* end of loop */
        leaq 64(%r10),%r10
        leaq 64(%r11),%r11
        addq $8,%rax
	cmpq %r9,%rax
	jb fftsize2loop

//    //size 4
//    //[1  0  1  0]
//    //[0  1  0 -i]
//    //[1  0 -1  0]
//    //[0  1  0  i]
//    // r0 + r2    i0 + i2
//    // r1 + i3    i1 - r3
//    // r0 - r2    i0 - i2
//    // r1 - i3    i1 + r3
//    {
//	for (int32_t block=0; block<ns4; block+=4) {
//	    double* re = pre+block;
//	    double* im = pim+block;
//	    tmp0[0]=re[0];
//	    tmp0[1]=re[1];
//	    tmp0[2]=re[0];
//	    tmp0[3]=re[1];
//	    tmp1[0]=re[2];
//	    tmp1[1]=im[3];
//	    tmp1[2]=-re[2];
//	    tmp1[3]=-im[3];
//	    tmp2[0]=im[0];
//	    tmp2[1]=im[1];
//	    tmp2[2]=im[0];
//	    tmp2[3]=im[1];
//	    tmp3[0]=im[2];
//	    tmp3[1]=-re[3];
//	    tmp3[2]=-im[2];
//	    tmp3[3]=re[3];
//	    add4(re,tmp0,tmp1);
//	    add4(im,tmp2,tmp3);
//	}
//    }
	movq	$0, %rax
	movq	%rdi,%r10
	movq	%rsi,%r11
fftsize4loop:
	vmovapd (%r10),%zmm0 /* r0 r1 r2 r3 */
	vmovapd (%r11),%zmm1 /* i0 i1 i2 i3 */
	// vperm2f128 $32,%zmm0,%zmm0,%zmm4 /* r0 r1 r0 r1 */
	// vperm2f128 $32,%zmm1,%zmm1,%zmm5 /* i0 i1 i0 i1 */
	vpermpd $68, %zmm0,%zmm4
	vpermpd $68, %zmm1,%zmm5
	# vperm2f128 $49,%zmm0,%zmm0,%zmm6 /* r2 r3 r2 r3 */
	# vperm2f128 $49,%zmm1,%zmm1,%zmm7 /* i2 i3 i2 i3 */
	# vshufpd $10,%zmm7,%zmm6,%zmm8    /* r2 i3 r2 i3 */
	# vshufpd $10,%zmm6,%zmm7,%zmm9    /* i2 r3 i2 r3 */
	vmovapd %zmm11, %zmm8
	vpermi2pd %zmm1, %zmm0, %zmm8
	vmovapd %zmm11, %zmm9
	vpermi2pd %zmm0, %zmm1, %zmm9

	vfmadd231pd %zmm8,%zmm13,%zmm4   /* (r0 r1 r0 r1) + (r2 i3 -r2 -i3) */
	vfmadd231pd %zmm9,%zmm14,%zmm5   /* (i0 i1 i0 i1) + (i2 -r3 -i2 r3) */
	vmovapd %zmm4,(%r10)
  vmovapd %zmm5,(%r11)
        /* end of loop */
        leaq 64(%r10),%r10
        leaq 64(%r11),%r11
        addq $8,%rax
	cmpq %r9,%rax
	jb fftsize4loop


//    
//    //general loop
//    const double* cur_tt = trig_tables;
//    for (int32_t halfnn=4; halfnn<ns4; halfnn*=2) {
//	int32_t nn = 2*halfnn;
//	for (int32_t block=0; block<ns4; block+=nn) {
//	    for (int32_t off=0; off<halfnn; off+=4) {
//		double* re0 = pre + block + off;
//		double* im0 = pim + block + off;
//		double* re1 = pre + block + halfnn + off;
//		double* im1 = pim + block + halfnn + off;
//		const double* tcs = cur_tt+2*off;
//		const double* tsn = tcs+4;
//		dotp4(tmp0,re1,tcs); // re*cos
//		dotp4(tmp1,re1,tsn); // re*sin
//		dotp4(tmp2,im1,tcs); // im*cos
//		dotp4(tmp3,im1,tsn); // im*sin
//		sub4(tmp0,tmp0,tmp3); // re2
//		add4(tmp1,tmp1,tmp2); // im2
//		add4(tmp2,re0,tmp0); // re + re
//		add4(tmp3,im0,tmp1); // im + im
//		sub4(tmp0,re0,tmp0); // re - re
//		sub4(tmp1,im0,tmp1); // im - im
//		copy4(re0,tmp2);
//		copy4(im0,tmp3);
//		copy4(re1,tmp0);
//		copy4(im1,tmp1);
//	    }
//	}
//	cur_tt += nn;
//    }

# size 8 pass: replace halfnn=4 YMM loop with ZMM, processing 2 blocks/iter
# Trig for halfnn=4: cos(-2pi*k/8) = [1, 1/sqrt2, 0, -1/sqrt2] for k=0..3
#                    sin(-2pi*k/8) = [0, -1/sqrt2, -1, -1/sqrt2]
# Load hardcoded ZMM constants (two repetitions for both blocks)
vmovapd fftsize8cos(%rip), %zmm8
vmovapd fftsize8sin(%rip), %zmm9
movq $0, %rax          /* rax: block index in doubles (step 16 = 2 blocks) */
fftsize8loop:
	# Load block B: pre[rax..rax+7] = [re0_B, re1_B]
	#              pim[rax..rax+7] = [im0_B, im1_B]
	# Load block B+8: pre[rax+8..rax+15] = [re0_{B+8}, re1_{B+8}]
	vmovapd    (%rdi,%rax,8), %zmm0  /* [re0_B, re1_B] */
	vmovapd 64(%rdi,%rax,8), %zmm1  /* [re0_{B+8}, re1_{B+8}] */
	vmovapd    (%rsi,%rax,8), %zmm2  /* [im0_B, im1_B] */
	vmovapd 64(%rsi,%rax,8), %zmm3  /* [im0_{B+8}, im1_{B+8}] */
	# Interleave: gather re0/re1/im0/im1 across both blocks
	vshuff64x2 $0x44, %zmm1, %zmm0, %zmm4  /* [re0_B, re0_{B+8}] */
	vshuff64x2 $0xEE, %zmm1, %zmm0, %zmm5  /* [re1_B, re1_{B+8}] */
	vshuff64x2 $0x44, %zmm3, %zmm2, %zmm6  /* [im0_B, im0_{B+8}] */
	vshuff64x2 $0xEE, %zmm3, %zmm2, %zmm7  /* [im1_B, im1_{B+8}] */
	# Butterfly: re2 = re1*cos - im1*sin,  im2 = re1*sin + im1*cos
	vmulpd       %zmm5, %zmm8, %zmm10       /* re1*cos */
	vfnmadd231pd %zmm7, %zmm9, %zmm10       /* re2 = re1*cos - im1*sin */
	vmulpd       %zmm5, %zmm9, %zmm11       /* re1*sin */
	vfmadd231pd  %zmm7, %zmm8, %zmm11       /* im2 = re1*sin + im1*cos */
	# Output
	vsubpd %zmm10, %zmm4, %zmm5  /* new_re1 = re0 - re2 */
	vsubpd %zmm11, %zmm6, %zmm7  /* new_im1 = im0 - im2 */
	vaddpd %zmm10, %zmm4, %zmm4  /* new_re0 = re0 + re2 */
	vaddpd %zmm11, %zmm6, %zmm6  /* new_im0 = im0 + im2 */
	# Deinterleave: recombine [new_re0_B, new_re1_B] and [new_re0_{B+8}, new_re1_{B+8}]
	vshuff64x2 $0x44, %zmm5, %zmm4, %zmm0  /* [new_re0_B, new_re1_B] */
	vshuff64x2 $0xEE, %zmm5, %zmm4, %zmm1  /* [new_re0_{B+8}, new_re1_{B+8}] */
	vshuff64x2 $0x44, %zmm7, %zmm6, %zmm2  /* [new_im0_B, new_im1_B] */
	vshuff64x2 $0xEE, %zmm7, %zmm6, %zmm3  /* [new_im0_{B+8}, new_im1_{B+8}] */
	# Store
	vmovapd %zmm0,    (%rdi,%rax,8)
	vmovapd %zmm1, 64(%rdi,%rax,8)
	vmovapd %zmm2,    (%rsi,%rax,8)
	vmovapd %zmm3, 64(%rsi,%rax,8)
	addq $16, %rax
	cmpq %r9, %rax
	jb fftsize8loop
# ── Option B: fused {halfnn=8, halfnn=16} pass ──────────────────────────────
# Trig table layout relative to %r8 (base of runtime trig table):
#   [+0  ..+63 ] halfnn=4 entry (unused by fftsize8loop; use hardcoded consts)
#   [+64 ..+191] halfnn=8 entry: W8_cos[0..7] at +64, W8_sin[0..7] at +128
#   [+192..+447] halfnn=16 entry:
#                  W16[off=0]_cos at +192, W16[off=0]_sin at +256
#                  W16[off=8]_cos at +320, W16[off=8]_sin at +384
#   [+448..    ] halfnn=32+ entries → start of general loop
#
# Super-block of 32 doubles [A=base+0..7, B=base+8..15, C=base+16..23, D=base+24..31]:
#   Step 1: halfnn=8 butterfly A<->B (trig W8)  and C<->D (same W8)
#   Step 2: halfnn=16 butterfly A'<->C' (trig W16[off=0])  and B'<->D' (trig W16[off=8])
#
# Register map: zmm0-3=re data, zmm8-11=im data, zmm4-5=tmp,
#               zmm12-13=W8 cos/sin, zmm14-17=W16[0]/W16[8] cos/sin
	cmpq	$32,%r9
	jb	fft_fallback_to_8
	vmovapd  64(%r8),%zmm12		/* W8_cos */
	vmovapd 128(%r8),%zmm13		/* W8_sin */
	vmovapd 192(%r8),%zmm14		/* W16[off=0]_cos */
	vmovapd 256(%r8),%zmm15		/* W16[off=0]_sin */
	vmovapd 320(%r8),%zmm16		/* W16[off=8]_cos */
	vmovapd 384(%r8),%zmm17		/* W16[off=8]_sin */
	movq	$0,%rbx			/* rbx: super-block base (doubles) */
.p2align 4
fft_r4_8_16_loop:
	vmovapd    (%rdi,%rbx,8),%zmm0		/* A_re */
	vmovapd  64(%rdi,%rbx,8),%zmm1		/* B_re */
	vmovapd 128(%rdi,%rbx,8),%zmm2		/* C_re */
	vmovapd 192(%rdi,%rbx,8),%zmm3		/* D_re */
	vmovapd    (%rsi,%rbx,8),%zmm8		/* A_im */
	vmovapd  64(%rsi,%rbx,8),%zmm9		/* B_im */
	vmovapd 128(%rsi,%rbx,8),%zmm10	/* C_im */
	vmovapd 192(%rsi,%rbx,8),%zmm11	/* D_im */
	# halfnn=8 butterfly: A vs B
	vmulpd	%zmm1,%zmm12,%zmm4		/* B_re*W8_cos */
	vmulpd	%zmm1,%zmm13,%zmm5		/* B_re*W8_sin */
	vfnmadd231pd %zmm9,%zmm13,%zmm4	/* re2 = B_re*cos - B_im*sin */
	vfmadd231pd  %zmm9,%zmm12,%zmm5	/* im2 = B_re*sin + B_im*cos */
	vsubpd	%zmm4,%zmm0,%zmm1		/* new_B_re = A_re - re2 */
	vsubpd	%zmm5,%zmm8,%zmm9		/* new_B_im = A_im - im2 */
	vaddpd	%zmm4,%zmm0,%zmm0		/* new_A_re = A_re + re2 */
	vaddpd	%zmm5,%zmm8,%zmm8		/* new_A_im = A_im + im2 */
	# halfnn=8 butterfly: C vs D
	vmulpd	%zmm3,%zmm12,%zmm4		/* D_re*W8_cos */
	vmulpd	%zmm3,%zmm13,%zmm5		/* D_re*W8_sin */
	vfnmadd231pd %zmm11,%zmm13,%zmm4	/* re2 = D_re*cos - D_im*sin */
	vfmadd231pd  %zmm11,%zmm12,%zmm5	/* im2 = D_re*sin + D_im*cos */
	vsubpd	%zmm4,%zmm2,%zmm3		/* new_D_re = C_re - re2 */
	vsubpd	%zmm5,%zmm10,%zmm11		/* new_D_im = C_im - im2 */
	vaddpd	%zmm4,%zmm2,%zmm2		/* new_C_re = C_re + re2 */
	vaddpd	%zmm5,%zmm10,%zmm10		/* new_C_im = C_im + im2 */
	# halfnn=16 butterfly: A' vs C'
	vmulpd	%zmm2,%zmm14,%zmm4		/* C_re*W16_0_cos */
	vmulpd	%zmm2,%zmm15,%zmm5		/* C_re*W16_0_sin */
	vfnmadd231pd %zmm10,%zmm15,%zmm4	/* re2 = C_re*cos - C_im*sin */
	vfmadd231pd  %zmm10,%zmm14,%zmm5	/* im2 = C_re*sin + C_im*cos */
	vsubpd	%zmm4,%zmm0,%zmm2		/* new_C_re = A_re - re2 */
	vsubpd	%zmm5,%zmm8,%zmm10		/* new_C_im = A_im - im2 */
	vaddpd	%zmm4,%zmm0,%zmm0		/* new_A_re = A_re + re2 */
	vaddpd	%zmm5,%zmm8,%zmm8		/* new_A_im = A_im + im2 */
	# halfnn=16 butterfly: B' vs D'
	vmulpd	%zmm3,%zmm16,%zmm4		/* D_re*W16_8_cos */
	vmulpd	%zmm3,%zmm17,%zmm5		/* D_re*W16_8_sin */
	vfnmadd231pd %zmm11,%zmm17,%zmm4	/* re2 = D_re*cos - D_im*sin */
	vfmadd231pd  %zmm11,%zmm16,%zmm5	/* im2 = D_re*sin + D_im*cos */
	vsubpd	%zmm4,%zmm1,%zmm3		/* new_D_re = B_re - re2 */
	vsubpd	%zmm5,%zmm9,%zmm11		/* new_D_im = B_im - im2 */
	vaddpd	%zmm4,%zmm1,%zmm1		/* new_B_re = B_re + re2 */
	vaddpd	%zmm5,%zmm9,%zmm9		/* new_B_im = B_im + im2 */
	vmovapd %zmm0,   (%rdi,%rbx,8)
	vmovapd %zmm1, 64(%rdi,%rbx,8)
	vmovapd %zmm2,128(%rdi,%rbx,8)
	vmovapd %zmm3,192(%rdi,%rbx,8)
	vmovapd %zmm8,   (%rsi,%rbx,8)
	vmovapd %zmm9, 64(%rsi,%rbx,8)
	vmovapd %zmm10,128(%rsi,%rbx,8)
	vmovapd %zmm11,192(%rsi,%rbx,8)
	addq	$32,%rbx
	cmpq	%r9,%rbx
	jb	fft_r4_8_16_loop
	# General loop for halfnn=32+ (skip halfnn=4,8,16 trig entries: 64+128+256=448 bytes)
	leaq	448(%r8),%rdx
	movq	$32,%rax
	cmpq	%r9,%rax
	jb	ffthalfnnloop
	jmp	fftbeforefinal
fft_fallback_to_8:
	leaq	64(%r8),%rdx		/* skip halfnn=4 entry; halfnn=8 at r8+64 */
	movq	$8,%rax
ffthalfnnloop:
	movq $0,%rbx /* rbx: block */
fftblockloop:
	leaq (%rdi,%rbx,8),%r10 /* re0 pointer */
	leaq (%rsi,%rbx,8),%r11 /* im0 pointer */
	leaq (%r10,%rax,8),%r12 /* re1 pointer */
	leaq (%r11,%rax,8),%r13 /* im1 pointer */
	movq %rdx,%r14          /* tcs pointer */
	movq $0,%rcx /* rcx: off */
fftoffloop:
	vmovapd (%r10),%zmm0 /* re0 */
	vmovapd (%r11),%zmm1 /* im0 */
	vmovapd (%r12),%zmm2 /* re1 */
	vmovapd (%r13),%zmm3 /* im1 */
	vmovapd (%r14),%zmm4 /* cos */
	vmovapd 64(%r14),%zmm5 /* sin */
	vmulpd	%zmm2,%zmm4,%zmm6 /* re1.cos */
	vmulpd	%zmm2,%zmm5,%zmm7 /* re1.sin */
        vfnmadd231pd %zmm3,%zmm5,%zmm6 /* re2 = re1.cos - im1.sin */
        vfmadd231pd %zmm3,%zmm4,%zmm7  /* im2 = re1.sin + im1.cos */
	vsubpd	%zmm6,%zmm0,%zmm2 /* re0 - re2 */
	vsubpd	%zmm7,%zmm1,%zmm3 /* im0 - im2 */
	vaddpd	%zmm6,%zmm0,%zmm0 /* re0 + re2 */
	vaddpd	%zmm7,%zmm1,%zmm1 /* im0 + im2 */
	vmovapd %zmm0,(%r10)
	vmovapd %zmm1,(%r11)
	vmovapd %zmm2,(%r12)
	vmovapd %zmm3,(%r13)
        /* end of off loop */
    	leaq 	64(%r10),%r10
    	leaq	64(%r11),%r11
    	leaq 	64(%r12),%r12
    	leaq 	64(%r13),%r13
	leaq 	128(%r14),%r14
	addq 	$8,%rcx
	cmpq	%rax,%rcx
	jb 	fftoffloop
	/* end of block loop */
	leaq	(%rbx,%rax,2),%rbx
	cmpq	%r9,%rbx
	jb 	fftblockloop
	/* end of halfnn loop */
	shlq	$1,%rax
	leaq	(%rdx,%rax,8),%rdx
	cmpq	%r9,%rax
	jb ffthalfnnloop
fftbeforefinal:

//    //multiply by omb^j
//    for (int32_t j=0; j<ns4; j+=4) {
//	const double* r0 = cur_tt+2*j;
//	const double* r1 = r0+4;
//	//(re*cos-im*sin) + i (im*cos+re*sin)
//	double* d0 = pre+j;
//	double* d1 = pim+j;
//	dotp4(tmp0,d0,r0); //re*cos
//	dotp4(tmp1,d1,r0); //im*cos
//	dotp4(tmp2,d0,r1); //re*sin
//	dotp4(tmp3,d1,r1); //im*sin
//	sub4(d0,tmp0,tmp3);
//	add4(d1,tmp1,tmp2);
//    }
//}
	/* cur_tt is at rdx */
	movq $0,%rax /* j */
	movq %rdi,%r10
	movq %rsi,%r11
fftfinalloop:
	vmovapd	(%r10),%zmm0 /* re */
	vmovapd	(%r11),%zmm1 /* im */
	vmovapd (%rdx),%zmm2 /* cos */
	vmovapd 64(%rdx),%zmm3 /* sin */
	vmulpd %zmm0,%zmm2,%zmm4        /* re*cos */
	vmulpd %zmm0,%zmm3,%zmm5        /* re*sin */
	vfnmadd231pd %zmm1,%zmm3,%zmm4  /* re*cos - im*sin */
	vfmadd231pd %zmm1,%zmm2,%zmm5   /* re*sin + im*cos */
	vmovapd %zmm4,%zmm0
	vmovapd %zmm5,%zmm1
	vmovapd %zmm0,(%r10)
	vmovapd %zmm1,(%r11)
    	/* end of final loop */
    	leaq	64(%r10),%r10
    	leaq	64(%r11),%r11
    	leaq	128(%rdx),%rdx
	addq	$8,%rax
	cmpq	%r9,%rax
	jb fftfinalloop

	/* Restore registers */
fftend:
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
size4negation0: .double +1.0, +1.0, +1.0, -1.0, +1.0, +1.0, +1.0, -1.0 /* zmm15 */
size4negation1: .double +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0, +1.0 /* zmm14 */
size4negation2: .double +1.0, +1.0, -1.0, -1.0, +1.0, +1.0, -1.0, -1.0 /* zmm13 */
size4negation3: .double +1.0, -1.0, +1.0, -1.0, +1.0, -1.0, +1.0, -1.0 /* zmm12 */
permutex1: .quad 2, 8+3, 2, 8+3, 6, 8+7, 6, 8+7 /* zmm11 */
/* FFT size-8 pass trig constants: cos/sin(-2pi*k/8) for k=0..3, repeated twice for ZMM */
.balign 64
fftsize8cos: .double +1.0, +0.7071067811865476, +0.0, -0.7071067811865476, +1.0, +0.7071067811865476, +0.0, -0.7071067811865476
fftsize8sin: .double +0.0, -0.7071067811865476, -1.0, -0.7071067811865476, +0.0, -0.7071067811865476, -1.0, -0.7071067811865476

#if !__APPLE__
	.size	fft, .-fft
#endif
