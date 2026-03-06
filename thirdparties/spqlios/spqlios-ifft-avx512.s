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
	cmpq	$256, %r10
	je	ifft_fused_twist_d_setup	/* N=1024: skip firstloop, fuse twist with Option D */
	movq    $0, %rcx  /* Loop counter: Range [0, r10), step size 4 */
	movq	%r8, %r11 /* r11 is the trig table pointer, step size 64 */
.p2align 5
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

# ── Option D: fused {halfnn=128, halfnn=64} pass ────────────────────────────
# DIF order: process halfnn=128 first (large), then halfnn=64 (small)
# 4-column groups: (i, i+64, i+128, i+192) for i=0..63 step 8
# Each super-block covers 256 elements. Trig repeats across super-blocks.
# Trig: W128[i] + W128[i+64] (4 ZMMs) + W64[i] (2 ZMMs) = 6 trig ZMMs
# Data: 4 re + 4 im = 8 ZMMs. Temp: zmm30-31. Total: 16 ZMMs.
	# Advance r8 for halfnn=128 (nn=256, r12=256)
	leaq (%r8,%r12,8),%r8
	leaq (%r8,%r12,8),%r8		/* r8 = W128 trig base */
	movq %r8,%rcx			/* rcx = W128 lo trig base */
	leaq 1024(%rcx),%r14		/* r14 = W128 hi trig base */
	# Advance r8 for halfnn=64 (nn=128)
	movq $128,%r13
	leaq (%r8,%r13,8),%r8
	leaq (%r8,%r13,8),%r8		/* r8 = W64 trig base */
	movq %r8,%rdx			/* rdx = W64 trig base */
	# Save trig bases for super-block reset
	# r9 is free in IFFT (not used as ns4 constant here; r10 is ns4)
	movq %rcx,%r9			/* r9 = saved W128 lo base */
	movq %r14,%r13			/* r13 = saved W128 hi base */
	movq %rdx,%r11			/* r11 = saved W64 base */
	movq $0,%rbx			/* rbx = data byte offset */
	leaq (,%r10,8),%rax		/* rax = ns4 * 8 = total data bytes */
.p2align 4
ifft_d128_64_outer:
	movq %r9,%rcx			/* reset W128 lo trig */
	movq %r13,%r14			/* reset W128 hi trig */
	movq %r11,%rdx			/* reset W64 trig */
	leaq 512(%rbx),%r12		/* r12 = end of this super-block */
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
	# Advance inner pointers
	addq	$64,%rbx		/* data: +8 doubles */
	leaq	128(%rcx),%rcx		/* W128 lo: +128 bytes */
	leaq	128(%r14),%r14		/* W128 hi: +128 bytes */
	leaq	128(%rdx),%rdx		/* W64: +128 bytes */
	cmpq	%r12,%rbx		/* end of this super-block? */
	jb	ifft_d128_64_loop
	# Advance to next super-block: skip from sb+512 to sb+2048
	leaq	1536(%rbx),%rbx
	cmpq	%rax,%rbx		/* done with all super-blocks? */
	jb	ifft_d128_64_outer
	# r8 is at W64 trig base (correct for Option C to advance from)
	movq	$64,%r12		/* r12 = 64 (nn for Option C) */
	jmp	ifft_option_c_entry

ifft_fused_twist_d_setup:
	# Fused twist + Option D for N=1024 (ns4=256)
	# Applies firstloop twist AND halfnn=128+64 DIF butterflies in one pass.
	# Eliminates the separate firstloop and Option D store/load round-trips.
	movq	%r8,%rax		/* rax = twist trig base (original r8) */
	movq	%r10,%r12		/* r12 = ns4 = 256 */
	# Advance r8 for halfnn=128 (nn=256)
	leaq	(%r8,%r12,8),%r8
	leaq	(%r8,%r12,8),%r8	/* r8 = W128 trig base */
	movq	%r8,%rcx		/* rcx = W128 trig pointer */
	leaq	1024(%rcx),%r14		/* r14 = W128[i+64] trig pointer */
	# Advance r8 for halfnn=64 (nn=128)
	movq	$128,%r13
	leaq	(%r8,%r13,8),%r8
	leaq	(%r8,%r13,8),%r8	/* r8 = W64 trig base */
	movq	%r8,%rdx		/* rdx = W64 trig pointer */
	movq	$0,%rbx
.p2align 4
ifft_fused_twist_d_loop:
	# Load 4 columns of data
	vmovapd	(%rdi,%rbx),%zmm0
	vmovapd	512(%rdi,%rbx),%zmm1
	vmovapd	1024(%rdi,%rbx),%zmm2
	vmovapd	1536(%rdi,%rbx),%zmm3
	vmovapd	(%rsi,%rbx),%zmm4
	vmovapd	512(%rsi,%rbx),%zmm5
	vmovapd	1024(%rsi,%rbx),%zmm6
	vmovapd	1536(%rsi,%rbx),%zmm7
	# ── Apply twist to each column in-register ──
	# Column 0 (zmm0=re, zmm4=im): twist at rax
	vmovapd	(%rax),%zmm8
	vmovapd	64(%rax),%zmm9
	vmulpd	%zmm0,%zmm9,%zmm10
	vmulpd	%zmm0,%zmm8,%zmm0
	vfnmadd231pd %zmm4,%zmm9,%zmm0
	vfmadd231pd  %zmm4,%zmm8,%zmm10
	vmovapd	%zmm10,%zmm4
	# Column 1 (zmm1=re, zmm5=im): twist at rax+1024
	vmovapd	1024(%rax),%zmm8
	vmovapd	1088(%rax),%zmm9
	vmulpd	%zmm1,%zmm9,%zmm10
	vmulpd	%zmm1,%zmm8,%zmm1
	vfnmadd231pd %zmm5,%zmm9,%zmm1
	vfmadd231pd  %zmm5,%zmm8,%zmm10
	vmovapd	%zmm10,%zmm5
	# Column 2 (zmm2=re, zmm6=im): twist at rax+2048
	vmovapd	2048(%rax),%zmm8
	vmovapd	2112(%rax),%zmm9
	vmulpd	%zmm2,%zmm9,%zmm10
	vmulpd	%zmm2,%zmm8,%zmm2
	vfnmadd231pd %zmm6,%zmm9,%zmm2
	vfmadd231pd  %zmm6,%zmm8,%zmm10
	vmovapd	%zmm10,%zmm6
	# Column 3 (zmm3=re, zmm7=im): twist at rax+3072
	vmovapd	3072(%rax),%zmm8
	vmovapd	3136(%rax),%zmm9
	vmulpd	%zmm3,%zmm9,%zmm10
	vmulpd	%zmm3,%zmm8,%zmm3
	vfnmadd231pd %zmm7,%zmm9,%zmm3
	vfmadd231pd  %zmm7,%zmm8,%zmm10
	vmovapd	%zmm10,%zmm7
	# ── IFFT DIF butterflies (identical to regular Option D) ──
	# halfnn=128 DIF: zmm0<->zmm2 with W128[i]
	vmovapd	(%rcx),%zmm8
	vmovapd	64(%rcx),%zmm9
	vsubpd	%zmm2,%zmm0,%zmm30
	vsubpd	%zmm6,%zmm4,%zmm31
	vaddpd	%zmm2,%zmm0,%zmm0
	vaddpd	%zmm6,%zmm4,%zmm4
	vmulpd	%zmm30,%zmm8,%zmm2
	vfnmadd231pd %zmm31,%zmm9,%zmm2
	vmulpd	%zmm30,%zmm9,%zmm6
	vfmadd231pd  %zmm31,%zmm8,%zmm6
	# halfnn=128 DIF: zmm1<->zmm3 with W128[i+64]
	vmovapd	(%r14),%zmm8
	vmovapd	64(%r14),%zmm9
	vsubpd	%zmm3,%zmm1,%zmm30
	vsubpd	%zmm7,%zmm5,%zmm31
	vaddpd	%zmm3,%zmm1,%zmm1
	vaddpd	%zmm7,%zmm5,%zmm5
	vmulpd	%zmm30,%zmm8,%zmm3
	vfnmadd231pd %zmm31,%zmm9,%zmm3
	vmulpd	%zmm30,%zmm9,%zmm7
	vfmadd231pd  %zmm31,%zmm8,%zmm7
	# halfnn=64 DIF: zmm0<->zmm1 with W64[i]
	vmovapd	(%rdx),%zmm8
	vmovapd	64(%rdx),%zmm9
	vsubpd	%zmm1,%zmm0,%zmm30
	vsubpd	%zmm5,%zmm4,%zmm31
	vaddpd	%zmm1,%zmm0,%zmm0
	vaddpd	%zmm5,%zmm4,%zmm4
	vmulpd	%zmm30,%zmm8,%zmm1
	vfnmadd231pd %zmm31,%zmm9,%zmm1
	vmulpd	%zmm30,%zmm9,%zmm5
	vfmadd231pd  %zmm31,%zmm8,%zmm5
	# halfnn=64 DIF: zmm2<->zmm3 with W64[i] (same trig)
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
	# Advance
	addq	$64,%rbx
	leaq	128(%rcx),%rcx
	leaq	128(%r14),%r14
	leaq	128(%rdx),%rdx
	leaq	128(%rax),%rax
	cmpq	$512,%rbx
	jb	ifft_fused_twist_d_loop
	# r8 is already at W64 trig base (set up before loop)
	movq	$64,%r12
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
# ── Fused size-8 + size-4 + size-2 pass ──────────────────────────────────────
# Replaces separate ifftsize8loop, size4loop, size2loop
# Processes 16 doubles (2 blocks of 8) per iteration, doing all 3 stages in-register.
# DIF order: size-8 (halfnn=4) → size-4 → size-2
#
# Register map:
#   zmm24 = ifftsize8cos (for size-8)
#   zmm25 = ifftsize8sin (for size-8)
#   zmm26 = permutex1 [0,1,0,9,4,5,4,13] (for size-4)
#   zmm27 = permutex2 [2,3,2,11,6,7,6,15] (for size-4)
#   zmm28 = size4negation0 [+1,+1,+1,-1,...] (for size-4 re vmulpd)
#   zmm29 = size4negation1 [+1,+1,-1,+1,...] (for size-4 re vfmadd)
#   zmm30 = size4negation2 [+1,+1,-1,-1,...] (for size-4 im vfmadd)
#   zmm31 = size4negation3 [+1,-1,+1,-1,...] (for size-2)
#   zmm0-3 = loaded data; zmm4-23 = temps
	vmovapd     ifftsize8cos(%rip), %zmm24
	vmovapd     ifftsize8sin(%rip), %zmm25
	vmovapd     permutex1(%rip), %zmm26
	vmovapd     permutex2(%rip), %zmm27
	vmovapd     size4negation0(%rip), %zmm28
	vmovapd     size4negation1(%rip), %zmm29
	vmovapd     size4negation2(%rip), %zmm30
	vmovapd     size4negation3(%rip), %zmm31
	movq $0, %r11
.p2align 5
ifft_fused_842_loop:
	# Load 2 blocks: re[0..7], re[8..15], im[0..7], im[8..15]
	vmovapd    (%rdi,%r11,8), %zmm0
	vmovapd 64(%rdi,%r11,8), %zmm1
	vmovapd    (%rsi,%r11,8), %zmm2
	vmovapd 64(%rsi,%r11,8), %zmm3
	# ── Size-8 (halfnn=4) IFFT DIF butterfly across 2 blocks ──
	vshuff64x2 $0x44, %zmm1, %zmm0, %zmm4
	vshuff64x2 $0xEE, %zmm1, %zmm0, %zmm5
	vshuff64x2 $0x44, %zmm3, %zmm2, %zmm6
	vshuff64x2 $0xEE, %zmm3, %zmm2, %zmm7
	vaddpd %zmm4, %zmm5, %zmm8
	vaddpd %zmm6, %zmm7, %zmm9
	vsubpd %zmm5, %zmm4, %zmm4
	vsubpd %zmm7, %zmm6, %zmm6
	vmulpd       %zmm4, %zmm24, %zmm5
	vfnmadd231pd %zmm6, %zmm25, %zmm5
	vmulpd       %zmm4, %zmm25, %zmm7
	vfmadd231pd  %zmm6, %zmm24, %zmm7
	vshuff64x2 $0x44, %zmm5, %zmm8, %zmm0
	vshuff64x2 $0xEE, %zmm5, %zmm8, %zmm1
	vshuff64x2 $0x44, %zmm7, %zmm9, %zmm2
	vshuff64x2 $0xEE, %zmm7, %zmm9, %zmm3
	# Data: re_lo=zmm0, re_hi=zmm1, im_lo=zmm2, im_hi=zmm3
	# ── Size-4 IFFT DIF butterfly (cross re/im within each ZMM pair) ──
	# Pair 1: (zmm0=re_lo, zmm2=im_lo) → (zmm4=new_re_lo, zmm6=new_im_lo)
	vmovapd %zmm26, %zmm4
	vpermi2pd %zmm2, %zmm0, %zmm4
	vmovapd %zmm27, %zmm5
	vpermi2pd %zmm2, %zmm0, %zmm5
	vmovapd %zmm26, %zmm6
	vpermi2pd %zmm0, %zmm2, %zmm6
	vmovapd %zmm27, %zmm7
	vpermi2pd %zmm0, %zmm2, %zmm7
	vmulpd %zmm4, %zmm28, %zmm4
	vfmadd231pd %zmm5, %zmm29, %zmm4
	vfmadd231pd %zmm7, %zmm30, %zmm6
	# Pair 2: (zmm1=re_hi, zmm3=im_hi) → (zmm8=new_re_hi, zmm10=new_im_hi)
	vmovapd %zmm26, %zmm8
	vpermi2pd %zmm3, %zmm1, %zmm8
	vmovapd %zmm27, %zmm9
	vpermi2pd %zmm3, %zmm1, %zmm9
	vmovapd %zmm26, %zmm10
	vpermi2pd %zmm1, %zmm3, %zmm10
	vmovapd %zmm27, %zmm11
	vpermi2pd %zmm1, %zmm3, %zmm11
	vmulpd %zmm8, %zmm28, %zmm8
	vfmadd231pd %zmm9, %zmm29, %zmm8
	vfmadd231pd %zmm11, %zmm30, %zmm10
	# Data: re_lo=zmm4, re_hi=zmm8, im_lo=zmm6, im_hi=zmm10
	# ── Size-2 butterfly (within each ZMM independently) ──
	vshufpd $0x00, %zmm4, %zmm4, %zmm0
	vshufpd $0xFF, %zmm4, %zmm4, %zmm4
	vfmadd231pd %zmm4, %zmm31, %zmm0
	vshufpd $0x00, %zmm8, %zmm8, %zmm1
	vshufpd $0xFF, %zmm8, %zmm8, %zmm8
	vfmadd231pd %zmm8, %zmm31, %zmm1
	vshufpd $0x00, %zmm6, %zmm6, %zmm2
	vshufpd $0xFF, %zmm6, %zmm6, %zmm6
	vfmadd231pd %zmm6, %zmm31, %zmm2
	vshufpd $0x00, %zmm10, %zmm10, %zmm3
	vshufpd $0xFF, %zmm10, %zmm10, %zmm10
	vfmadd231pd %zmm10, %zmm31, %zmm3
	# Store
	vmovapd %zmm0,    (%rdi,%r11,8)
	vmovapd %zmm1, 64(%rdi,%r11,8)
	vmovapd %zmm2,    (%rsi,%r11,8)
	vmovapd %zmm3, 64(%rsi,%r11,8)
	addq $16, %r11
	cmpq %r10, %r11
	jb ifft_fused_842_loop


	/* Restore registers */
end:
	vzeroupper
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

