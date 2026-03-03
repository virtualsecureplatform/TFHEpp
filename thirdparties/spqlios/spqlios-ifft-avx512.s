	.file	"spqlios-ifft-avx.s"
#if !__APPLE__
	.section .note.GNU-stack,"",%progbits
#endif
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
	cmpq $256,%r12
	jbe ifft_post_nnloop	/* skip nnloop for N <= 1024 */
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
	cmpq $256,%r12
	ja nnloop		/* loop while r12 > 256 */

ifft_post_nnloop:
	# r12 = 256 for N >= 1024, or ns4 for N <= 1024
	cmpq $256,%r12
	jne ifft_d_small_n

# ── Option D: fused {halfnn=128, halfnn=64} pass for N >= 1024 ──────────────
# DIF order: process halfnn=128 first (large), then halfnn=64 (small)
# 4-column groups: (i, i+64, i+128, i+192) for i=0..63 step 8
# Trig: W128[i] + W128[i+64] (4 ZMMs) + W64[i] (2 ZMMs) = 6 trig ZMMs
# Data: 4 re + 4 im = 8 ZMMs. Temp: zmm30-31. Total: 16 ZMMs.
	# Advance r8 for halfnn=128 (nn=256, r12=256)
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8		/* r8 = W128 trig base */
	movq %r8,%rcx			/* rcx = W128 trig pointer */
	leaq 1024(%rcx),%r14		/* r14 = W128[i+64] trig pointer */
	# Advance r8 for halfnn=64 (nn=128)
	movq $128,%r13
	leaq (%r8,%r13,8),%r8
	leaq (%r8,%r13,8),%r8		/* r8 = W64 trig base */
	movq %r8,%rdx			/* rdx = W64 trig pointer */
	movq $0,%rbx			/* rbx = byte offset */
.p2align 4
ifft_d128_64_loop:
	# Load 4 columns of data
	vmovapd	(%rdi,%rbx),%zmm0		/* re[i] */
	vmovapd	512(%rdi,%rbx),%zmm1		/* re[i+64] */
	vmovapd	1024(%rdi,%rbx),%zmm2		/* re[i+128] */
	vmovapd	1536(%rdi,%rbx),%zmm3		/* re[i+192] */
	vmovapd	(%rsi,%rbx),%zmm4		/* im[i] */
	vmovapd	512(%rsi,%rbx),%zmm5		/* im[i+64] */
	vmovapd	1024(%rsi,%rbx),%zmm6		/* im[i+128] */
	vmovapd	1536(%rsi,%rbx),%zmm7		/* im[i+192] */
	# Load W128[i] trig
	vmovapd	(%rcx),%zmm8			/* W128[i] cos */
	vmovapd	64(%rcx),%zmm9			/* W128[i] sin */
	# halfnn=128 DIF: zmm0 <-> zmm2 with W128[i]
	vsubpd	%zmm2,%zmm0,%zmm30
	vsubpd	%zmm6,%zmm4,%zmm31
	vaddpd	%zmm2,%zmm0,%zmm0
	vaddpd	%zmm6,%zmm4,%zmm4
	vmulpd	%zmm30,%zmm8,%zmm2
	vfnmadd231pd %zmm31,%zmm9,%zmm2
	vmulpd	%zmm30,%zmm9,%zmm6
	vfmadd231pd  %zmm31,%zmm8,%zmm6
	# Load W128[i+64] trig
	vmovapd	(%r14),%zmm8			/* W128[i+64] cos */
	vmovapd	64(%r14),%zmm9			/* W128[i+64] sin */
	# halfnn=128 DIF: zmm1 <-> zmm3 with W128[i+64]
	vsubpd	%zmm3,%zmm1,%zmm30
	vsubpd	%zmm7,%zmm5,%zmm31
	vaddpd	%zmm3,%zmm1,%zmm1
	vaddpd	%zmm7,%zmm5,%zmm5
	vmulpd	%zmm30,%zmm8,%zmm3
	vfnmadd231pd %zmm31,%zmm9,%zmm3
	vmulpd	%zmm30,%zmm9,%zmm7
	vfmadd231pd  %zmm31,%zmm8,%zmm7
	# Load W64[i] trig (same for both halfnn=64 pairs)
	vmovapd	(%rdx),%zmm8			/* W64[i] cos */
	vmovapd	64(%rdx),%zmm9			/* W64[i] sin */
	# halfnn=64 DIF: zmm0 <-> zmm1 with W64[i]
	vsubpd	%zmm1,%zmm0,%zmm30
	vsubpd	%zmm5,%zmm4,%zmm31
	vaddpd	%zmm1,%zmm0,%zmm0
	vaddpd	%zmm5,%zmm4,%zmm4
	vmulpd	%zmm30,%zmm8,%zmm1
	vfnmadd231pd %zmm31,%zmm9,%zmm1
	vmulpd	%zmm30,%zmm9,%zmm5
	vfmadd231pd  %zmm31,%zmm8,%zmm5
	# halfnn=64 DIF: zmm2 <-> zmm3 with W64[i] (same trig)
	vsubpd	%zmm3,%zmm2,%zmm30
	vsubpd	%zmm7,%zmm6,%zmm31
	vaddpd	%zmm3,%zmm2,%zmm2
	vaddpd	%zmm7,%zmm6,%zmm6
	vmulpd	%zmm30,%zmm8,%zmm3
	vfnmadd231pd %zmm31,%zmm9,%zmm3
	vmulpd	%zmm30,%zmm9,%zmm7
	vfmadd231pd  %zmm31,%zmm8,%zmm7
	# Store results
	vmovapd	%zmm0,(%rdi,%rbx)
	vmovapd	%zmm1,512(%rdi,%rbx)
	vmovapd	%zmm2,1024(%rdi,%rbx)
	vmovapd	%zmm3,1536(%rdi,%rbx)
	vmovapd	%zmm4,(%rsi,%rbx)
	vmovapd	%zmm5,512(%rsi,%rbx)
	vmovapd	%zmm6,1024(%rsi,%rbx)
	vmovapd	%zmm7,1536(%rsi,%rbx)
	# Advance pointers
	addq	$64,%rbx		/* data: +8 doubles */
	leaq	128(%rcx),%rcx		/* W128 lo: +128 bytes */
	leaq	128(%r14),%r14		/* W128 hi: +128 bytes */
	leaq	128(%rdx),%rdx		/* W64: +128 bytes */
	cmpq	$512,%rbx		/* done when rbx = 64*8 = 512 */
	jb	ifft_d128_64_loop
	# r8 is at W64 trig base (correct for Option C to advance from)
	movq	$64,%r12		/* r12 = 64 (nn for Option C) */
	jmp	ifft_option_c_entry

ifft_d_small_n:
	# N <= 512: process remaining stages individually until r12 < 128
	cmpq	$128,%r12
	jb	ifft_option_c_entry
ifft_d_small_loop:
	movq %r12,%r13
	shr  $1,%r13
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8
	movq $0,%r11
ifft_d_small_blockloop:
	leaq (%rdi,%r11,8),%rax
	leaq (%rsi,%r11,8),%rbx
	leaq (%rax,%r13,8),%rcx
	leaq (%rbx,%r13,8),%rdx
	movq $0,%r9
	movq %r8,%r14
ifft_d_small_offloop:
	vmovapd (%rax,%r9,8),%zmm0
	vmovapd (%rbx,%r9,8),%zmm1
	vmovapd (%rcx,%r9,8),%zmm2
	vmovapd (%rdx,%r9,8),%zmm3
	vaddpd %zmm0,%zmm2,%zmm4
	vaddpd %zmm1,%zmm3,%zmm5
	vsubpd %zmm2,%zmm0,%zmm6
	vsubpd %zmm3,%zmm1,%zmm7
	vmovapd %zmm4,(%rax,%r9,8)
	vmovapd %zmm5,(%rbx,%r9,8)
	vmovapd (%r14),%zmm8
	vmovapd 64(%r14),%zmm9
	vmulpd %zmm6,%zmm8,%zmm4
	vfnmadd231pd %zmm7,%zmm9,%zmm4
	vmulpd %zmm6,%zmm9,%zmm5
	vfmadd231pd %zmm7,%zmm8,%zmm5
	vmovapd %zmm4,(%rcx,%r9,8)
	vmovapd %zmm5,(%rdx,%r9,8)
	leaq 128(%r14),%r14
	addq $8,%r9
	cmpq %r13,%r9
	jb ifft_d_small_offloop
	addq %r12,%r11
	cmpq %r10,%r11
	jb ifft_d_small_blockloop
	movq %r13,%r12
	cmpq $128,%r12
	jae ifft_d_small_loop
ifft_option_c_entry:

# ── Option C: fused {halfnn=32, halfnn=16, halfnn=8} pass for N >= 512 ──────
# After nnloop/Option D (r12 = 64 for N >= 512):
#   r12=64: halfnn=32,16,8 unprocessed → fused pass
#   r12<64: fallback for smaller N
#
# IFFT DIF butterfly: new_re0=re0+re1, new_im0=im0+im1,
#   diff=(re0-re1, im0-im1); new_re1=diff_re*cos-diff_im*sin,
#                             new_im1=diff_re*sin+diff_im*cos
#
# Super-block of 64 doubles [A=+0,B=+8,C=+16,D=+24,E=+32,F=+40,G=+48,H=+56]
#   Step 1 (halfnn=32): A<->E (W32[0]), B<->F (W32[8]), C<->G (W32[16]), D<->H (W32[24])
#   Step 2 (halfnn=16): A'<->C' (W16[0]), B'<->D' (W16[8]), E'<->G' (W16[0]), F'<->H' (W16[8])
#   Step 3 (halfnn=8):  A''<->B'', C''<->D'', E''<->F'', G''<->H''  (all W8[0])
#
# Trig advances (from current r8 after nnloop exit with r12=64):
#   halfnn=32: advance r8 by 2*64*8=1024  (r8 → W32 base)
#   halfnn=16: advance r8 by 2*32*8=512   (r8 → W16 base)
#   halfnn=8:  advance r8 by 2*16*8=256   (r8 → W8 base)
#
# Register map: zmm0-7=re, zmm8-15=im, zmm16-17=W8, zmm18-21=W16, zmm22-29=W32, zmm30-31=tmp
	cmpq	$64,%r12
	jne	ifft_fallback_last
	# Advance r8 for halfnn=32 (nn=64, r12=64)
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8		/* r8 = W32 trig base */
	vmovapd   0(%r8),%zmm22		/* W32[0]_cos */
	vmovapd  64(%r8),%zmm23		/* W32[0]_sin */
	vmovapd 128(%r8),%zmm24		/* W32[8]_cos */
	vmovapd 192(%r8),%zmm25		/* W32[8]_sin */
	vmovapd 256(%r8),%zmm26		/* W32[16]_cos */
	vmovapd 320(%r8),%zmm27		/* W32[16]_sin */
	vmovapd 384(%r8),%zmm28		/* W32[24]_cos */
	vmovapd 448(%r8),%zmm29		/* W32[24]_sin */
	# Advance r8 for halfnn=16 (nn=32)
	movq	$32,%rcx
	leaq (%r8,%rcx,8),%r8
	leaq (%r8,%rcx,8),%r8		/* r8 = W16 trig base */
	vmovapd   0(%r8),%zmm18		/* W16[0]_cos */
	vmovapd  64(%r8),%zmm19		/* W16[0]_sin */
	vmovapd 128(%r8),%zmm20		/* W16[8]_cos */
	vmovapd 192(%r8),%zmm21		/* W16[8]_sin */
	# Advance r8 for halfnn=8 (nn=16)
	movq	$16,%rcx
	leaq (%r8,%rcx,8),%r8
	leaq (%r8,%rcx,8),%r8		/* r8 = W8 trig base */
	vmovapd   0(%r8),%zmm16		/* W8[0]_cos */
	vmovapd  64(%r8),%zmm17		/* W8[0]_sin */
	movq	$0,%r11			/* r11: super-block base (doubles) */
.p2align 4
ifft_r4_32_16_8_loop:
	vmovapd    (%rdi,%r11,8),%zmm0	/* A_re */
	vmovapd  64(%rdi,%r11,8),%zmm1	/* B_re */
	vmovapd 128(%rdi,%r11,8),%zmm2	/* C_re */
	vmovapd 192(%rdi,%r11,8),%zmm3	/* D_re */
	vmovapd 256(%rdi,%r11,8),%zmm4	/* E_re */
	vmovapd 320(%rdi,%r11,8),%zmm5	/* F_re */
	vmovapd 384(%rdi,%r11,8),%zmm6	/* G_re */
	vmovapd 448(%rdi,%r11,8),%zmm7	/* H_re */
	vmovapd    (%rsi,%r11,8),%zmm8	/* A_im */
	vmovapd  64(%rsi,%r11,8),%zmm9	/* B_im */
	vmovapd 128(%rsi,%r11,8),%zmm10	/* C_im */
	vmovapd 192(%rsi,%r11,8),%zmm11	/* D_im */
	vmovapd 256(%rsi,%r11,8),%zmm12	/* E_im */
	vmovapd 320(%rsi,%r11,8),%zmm13	/* F_im */
	vmovapd 384(%rsi,%r11,8),%zmm14	/* G_im */
	vmovapd 448(%rsi,%r11,8),%zmm15	/* H_im */
	# halfnn=32: A<->E (W32[0])
	vsubpd	%zmm4,%zmm0,%zmm30
	vsubpd	%zmm12,%zmm8,%zmm31
	vaddpd	%zmm4,%zmm0,%zmm0
	vaddpd	%zmm12,%zmm8,%zmm8
	vmulpd	%zmm30,%zmm22,%zmm4
	vfnmadd231pd %zmm31,%zmm23,%zmm4
	vmulpd	%zmm30,%zmm23,%zmm12
	vfmadd231pd  %zmm31,%zmm22,%zmm12
	# halfnn=32: B<->F (W32[8])
	vsubpd	%zmm5,%zmm1,%zmm30
	vsubpd	%zmm13,%zmm9,%zmm31
	vaddpd	%zmm5,%zmm1,%zmm1
	vaddpd	%zmm13,%zmm9,%zmm9
	vmulpd	%zmm30,%zmm24,%zmm5
	vfnmadd231pd %zmm31,%zmm25,%zmm5
	vmulpd	%zmm30,%zmm25,%zmm13
	vfmadd231pd  %zmm31,%zmm24,%zmm13
	# halfnn=32: C<->G (W32[16])
	vsubpd	%zmm6,%zmm2,%zmm30
	vsubpd	%zmm14,%zmm10,%zmm31
	vaddpd	%zmm6,%zmm2,%zmm2
	vaddpd	%zmm14,%zmm10,%zmm10
	vmulpd	%zmm30,%zmm26,%zmm6
	vfnmadd231pd %zmm31,%zmm27,%zmm6
	vmulpd	%zmm30,%zmm27,%zmm14
	vfmadd231pd  %zmm31,%zmm26,%zmm14
	# halfnn=32: D<->H (W32[24])
	vsubpd	%zmm7,%zmm3,%zmm30
	vsubpd	%zmm15,%zmm11,%zmm31
	vaddpd	%zmm7,%zmm3,%zmm3
	vaddpd	%zmm15,%zmm11,%zmm11
	vmulpd	%zmm30,%zmm28,%zmm7
	vfnmadd231pd %zmm31,%zmm29,%zmm7
	vmulpd	%zmm30,%zmm29,%zmm15
	vfmadd231pd  %zmm31,%zmm28,%zmm15
	# halfnn=16: A'<->C' (W16[0])
	vsubpd	%zmm2,%zmm0,%zmm30
	vsubpd	%zmm10,%zmm8,%zmm31
	vaddpd	%zmm2,%zmm0,%zmm0
	vaddpd	%zmm10,%zmm8,%zmm8
	vmulpd	%zmm30,%zmm18,%zmm2
	vfnmadd231pd %zmm31,%zmm19,%zmm2
	vmulpd	%zmm30,%zmm19,%zmm10
	vfmadd231pd  %zmm31,%zmm18,%zmm10
	# halfnn=16: B'<->D' (W16[8])
	vsubpd	%zmm3,%zmm1,%zmm30
	vsubpd	%zmm11,%zmm9,%zmm31
	vaddpd	%zmm3,%zmm1,%zmm1
	vaddpd	%zmm11,%zmm9,%zmm9
	vmulpd	%zmm30,%zmm20,%zmm3
	vfnmadd231pd %zmm31,%zmm21,%zmm3
	vmulpd	%zmm30,%zmm21,%zmm11
	vfmadd231pd  %zmm31,%zmm20,%zmm11
	# halfnn=16: E'<->G' (W16[0])
	vsubpd	%zmm6,%zmm4,%zmm30
	vsubpd	%zmm14,%zmm12,%zmm31
	vaddpd	%zmm6,%zmm4,%zmm4
	vaddpd	%zmm14,%zmm12,%zmm12
	vmulpd	%zmm30,%zmm18,%zmm6
	vfnmadd231pd %zmm31,%zmm19,%zmm6
	vmulpd	%zmm30,%zmm19,%zmm14
	vfmadd231pd  %zmm31,%zmm18,%zmm14
	# halfnn=16: F'<->H' (W16[8])
	vsubpd	%zmm7,%zmm5,%zmm30
	vsubpd	%zmm15,%zmm13,%zmm31
	vaddpd	%zmm7,%zmm5,%zmm5
	vaddpd	%zmm15,%zmm13,%zmm13
	vmulpd	%zmm30,%zmm20,%zmm7
	vfnmadd231pd %zmm31,%zmm21,%zmm7
	vmulpd	%zmm30,%zmm21,%zmm15
	vfmadd231pd  %zmm31,%zmm20,%zmm15
	# halfnn=8: A''<->B''
	vsubpd	%zmm1,%zmm0,%zmm30
	vsubpd	%zmm9,%zmm8,%zmm31
	vaddpd	%zmm1,%zmm0,%zmm0
	vaddpd	%zmm9,%zmm8,%zmm8
	vmulpd	%zmm30,%zmm16,%zmm1
	vfnmadd231pd %zmm31,%zmm17,%zmm1
	vmulpd	%zmm30,%zmm17,%zmm9
	vfmadd231pd  %zmm31,%zmm16,%zmm9
	# halfnn=8: C''<->D''
	vsubpd	%zmm3,%zmm2,%zmm30
	vsubpd	%zmm11,%zmm10,%zmm31
	vaddpd	%zmm3,%zmm2,%zmm2
	vaddpd	%zmm11,%zmm10,%zmm10
	vmulpd	%zmm30,%zmm16,%zmm3
	vfnmadd231pd %zmm31,%zmm17,%zmm3
	vmulpd	%zmm30,%zmm17,%zmm11
	vfmadd231pd  %zmm31,%zmm16,%zmm11
	# halfnn=8: E''<->F''
	vsubpd	%zmm5,%zmm4,%zmm30
	vsubpd	%zmm13,%zmm12,%zmm31
	vaddpd	%zmm5,%zmm4,%zmm4
	vaddpd	%zmm13,%zmm12,%zmm12
	vmulpd	%zmm30,%zmm16,%zmm5
	vfnmadd231pd %zmm31,%zmm17,%zmm5
	vmulpd	%zmm30,%zmm17,%zmm13
	vfmadd231pd  %zmm31,%zmm16,%zmm13
	# halfnn=8: G''<->H''
	vsubpd	%zmm7,%zmm6,%zmm30
	vsubpd	%zmm15,%zmm14,%zmm31
	vaddpd	%zmm7,%zmm6,%zmm6
	vaddpd	%zmm15,%zmm14,%zmm14
	vmulpd	%zmm30,%zmm16,%zmm7
	vfnmadd231pd %zmm31,%zmm17,%zmm7
	vmulpd	%zmm30,%zmm17,%zmm15
	vfmadd231pd  %zmm31,%zmm16,%zmm15
	vmovapd %zmm0,   (%rdi,%r11,8)
	vmovapd %zmm1, 64(%rdi,%r11,8)
	vmovapd %zmm2,128(%rdi,%r11,8)
	vmovapd %zmm3,192(%rdi,%r11,8)
	vmovapd %zmm4,256(%rdi,%r11,8)
	vmovapd %zmm5,320(%rdi,%r11,8)
	vmovapd %zmm6,384(%rdi,%r11,8)
	vmovapd %zmm7,448(%rdi,%r11,8)
	vmovapd %zmm8,   (%rsi,%r11,8)
	vmovapd %zmm9, 64(%rsi,%r11,8)
	vmovapd %zmm10,128(%rsi,%r11,8)
	vmovapd %zmm11,192(%rsi,%r11,8)
	vmovapd %zmm12,256(%rsi,%r11,8)
	vmovapd %zmm13,320(%rsi,%r11,8)
	vmovapd %zmm14,384(%rsi,%r11,8)
	vmovapd %zmm15,448(%rsi,%r11,8)
	addq	$64,%r11
	cmpq	%r10,%r11
	jb	ifft_r4_32_16_8_loop
	jmp	ifft_before_size8
ifft_fallback_last:
	# r12 != 64: run remaining iterations one at a time (small N / edge cases)
	# r12 ∈ {8,16,32} — loop until halfnn=4 (handled by ifftsize8loop)
ifft_fallback_iter:
	cmpq	$8,%r12
	jbe	ifft_before_size8
	movq %r12,%r13
	shr  $1,%r13
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8
	movq $0,%r11
ifft_fallback_blockloop:
	leaq (%rdi,%r11,8),%rax
	leaq (%rsi,%r11,8),%rbx
	leaq (%rax,%r13,8),%rcx
	leaq (%rbx,%r13,8),%rdx
	movq $0,%r9
	movq %r8,%r14
ifft_fallback_offloop:
	vmovapd (%rax,%r9,8),%zmm0
	vmovapd (%rbx,%r9,8),%zmm1
	vmovapd (%rcx,%r9,8),%zmm2
	vmovapd (%rdx,%r9,8),%zmm3
	vaddpd %zmm0,%zmm2,%zmm4
	vaddpd %zmm1,%zmm3,%zmm5
	vsubpd %zmm2,%zmm0,%zmm6
	vsubpd %zmm3,%zmm1,%zmm7
	vmovapd %zmm4,(%rax,%r9,8)
	vmovapd %zmm5,(%rbx,%r9,8)
	vmovapd (%r14),%zmm8
	vmovapd 64(%r14),%zmm9
	vmulpd %zmm6,%zmm8,%zmm4
	vfnmadd231pd %zmm7,%zmm9,%zmm4
	vmulpd %zmm6,%zmm9,%zmm5
	vfmadd231pd %zmm7,%zmm8,%zmm5
	vmovapd %zmm4,(%rcx,%r9,8)
	vmovapd %zmm5,(%rdx,%r9,8)
	leaq 128(%r14),%r14
	addq $8,%r9
	cmpq %r13,%r9
	jb ifft_fallback_offloop
	addq %r12,%r11
	cmpq %r10,%r11
	jb ifft_fallback_blockloop
	movq %r13,%r12
	jmp ifft_fallback_iter
ifft_before_size8:
	movq $0,%r11   /* r11 (block) */
# size 8 pass (IFFT last iteration): replace halfnn=4 YMM loop with ZMM, 2 blocks/iter
# Trig: conjugate of FFT: cos=[1, 1/sqrt2, 0, -1/sqrt2], sin=[0, +1/sqrt2, +1, +1/sqrt2]
vmovapd ifftsize8cos(%rip), %zmm8
vmovapd ifftsize8sin(%rip), %zmm9
movq $0, %r11  /* r11: block index in doubles (step 16 = 2 blocks) */
ifftsize8loop:
	# Load block B: are[r11..r11+7]=[re0_B,re1_B], block B+8: are[r11+8..r11+15]
	vmovapd    (%rdi,%r11,8), %zmm0  /* [re0_B, re1_B] */
	vmovapd 64(%rdi,%r11,8), %zmm1  /* [re0_{B+8}, re1_{B+8}] */
	vmovapd    (%rsi,%r11,8), %zmm2  /* [im0_B, im1_B] */
	vmovapd 64(%rsi,%r11,8), %zmm3  /* [im0_{B+8}, im1_{B+8}] */
	# Interleave: gather re0/re1/im0/im1 across both blocks
	vshuff64x2 $0x44, %zmm1, %zmm0, %zmm4  /* [re0_B, re0_{B+8}] */
	vshuff64x2 $0xEE, %zmm1, %zmm0, %zmm5  /* [re1_B, re1_{B+8}] */
	vshuff64x2 $0x44, %zmm3, %zmm2, %zmm6  /* [im0_B, im0_{B+8}] */
	vshuff64x2 $0xEE, %zmm3, %zmm2, %zmm7  /* [im1_B, im1_{B+8}] */
	# IFFT butterfly: new_re0=re0+re1, new_im0=im0+im1,
	#   new_re1=(re0-re1)*cos-(im0-im1)*sin, new_im1=(re0-re1)*sin+(im0-im1)*cos
	vaddpd %zmm4, %zmm5, %zmm10  /* new_re0 = re0 + re1 */
	vaddpd %zmm6, %zmm7, %zmm11  /* new_im0 = im0 + im1 */
	vsubpd %zmm5, %zmm4, %zmm4   /* diff_re = re0 - re1 */
	vsubpd %zmm7, %zmm6, %zmm6   /* diff_im = im0 - im1 */
	vmulpd       %zmm4, %zmm8, %zmm5       /* diff_re * cos */
	vfnmadd231pd %zmm6, %zmm9, %zmm5       /* new_re1 = diff_re*cos - diff_im*sin */
	vmulpd       %zmm4, %zmm9, %zmm7       /* diff_re * sin */
	vfmadd231pd  %zmm6, %zmm8, %zmm7       /* new_im1 = diff_re*sin + diff_im*cos */
	# Deinterleave: pack [new_re0_B, new_re1_B] and [new_re0_{B+8}, new_re1_{B+8}]
	vshuff64x2 $0x44, %zmm5, %zmm10, %zmm0  /* [new_re0_B, new_re1_B] */
	vshuff64x2 $0xEE, %zmm5, %zmm10, %zmm1  /* [new_re0_{B+8}, new_re1_{B+8}] */
	vshuff64x2 $0x44, %zmm7, %zmm11, %zmm2  /* [new_im0_B, new_im1_B] */
	vshuff64x2 $0xEE, %zmm7, %zmm11, %zmm3  /* [new_im0_{B+8}, new_im1_{B+8}] */
	# Store
	vmovapd %zmm0,    (%rdi,%r11,8)
	vmovapd %zmm1, 64(%rdi,%r11,8)
	vmovapd %zmm2,    (%rsi,%r11,8)
	vmovapd %zmm3, 64(%rsi,%r11,8)
	addq $16, %r11
	cmpq %r10, %r11
	jb ifftsize8loop


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
/* IFFT size-8 pass trig constants: cos/sin(+2pi*k/8) for k=0..3, repeated twice for ZMM */
.balign 64
ifftsize8cos: .double +1.0, +0.7071067811865476, +0.0, -0.7071067811865476, +1.0, +0.7071067811865476, +0.0, -0.7071067811865476
ifftsize8sin: .double +0.0, +0.7071067811865476, +1.0, +0.7071067811865476, +0.0, +0.7071067811865476, +1.0, +0.7071067811865476

#if !__APPLE__
	.size	ifft, .-ifft
#endif

