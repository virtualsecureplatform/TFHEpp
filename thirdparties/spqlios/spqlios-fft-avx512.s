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
# ── Fused size-2 + size-4 + size-8 pass ──────────────────────────────────────
# Replaces separate fftsize2loop, fftsize4loop, fftsize8loop
# Processes 16 doubles (2 blocks of 8) per iteration, doing all 3 stages in-register.
# DIT order: size-2 → size-4 → size-8 (halfnn=4)
#
# Register map:
#   zmm24 = size4negation3 [+1,-1,+1,-1,...] (for size-2)
#   zmm25 = size4negation2 [+1,+1,-1,-1,...] (for size-4 re)
#   zmm26 = size4negation1 [+1,-1,-1,+1,...] (for size-4 im)
#   zmm27 = permutex1 index (for size-4)
#   zmm28 = fftsize8cos (for size-8)
#   zmm29 = fftsize8sin (for size-8)
#   zmm0-3 = loaded data; zmm4-15 = temps
	vmovapd     size4negation3(%rip), %zmm24
	vmovapd     size4negation2(%rip), %zmm25
	vmovapd     size4negation1(%rip), %zmm26
	vmovapd     permutex1(%rip), %zmm27
	vmovapd     fftsize8cos(%rip), %zmm28
	vmovapd     fftsize8sin(%rip), %zmm29
	movq	$0,%rax
.p2align 5
fft_fused_248_loop:
	# Load 2 blocks: re[0..7], re[8..15], im[0..7], im[8..15]
	vmovapd    (%rdi,%rax,8), %zmm0
	vmovapd 64(%rdi,%rax,8), %zmm1
	vmovapd    (%rsi,%rax,8), %zmm2
	vmovapd 64(%rsi,%rax,8), %zmm3
	# ── Size-2 butterfly (within each ZMM independently) ──
	# [a0+a1, a0-a1, a2+a3, a2-a3, a4+a5, a4-a5, a6+a7, a6-a7]
	vshufpd $0x00, %zmm0, %zmm0, %zmm4
	vshufpd $0xFF, %zmm0, %zmm0, %zmm0
	vfmadd231pd %zmm0, %zmm24, %zmm4
	vshufpd $0x00, %zmm1, %zmm1, %zmm5
	vshufpd $0xFF, %zmm1, %zmm1, %zmm1
	vfmadd231pd %zmm1, %zmm24, %zmm5
	vshufpd $0x00, %zmm2, %zmm2, %zmm6
	vshufpd $0xFF, %zmm2, %zmm2, %zmm2
	vfmadd231pd %zmm2, %zmm24, %zmm6
	vshufpd $0x00, %zmm3, %zmm3, %zmm7
	vshufpd $0xFF, %zmm3, %zmm3, %zmm3
	vfmadd231pd %zmm3, %zmm24, %zmm7
	# Data: re_lo=zmm4, re_hi=zmm5, im_lo=zmm6, im_hi=zmm7
	# ── Size-4 butterfly (cross re/im within each ZMM pair) ──
	# Pair 1: (zmm4=re_lo, zmm6=im_lo) → (zmm0=new_re_lo, zmm2=new_im_lo)
	vpermpd $0x44, %zmm4, %zmm0
	vmovapd %zmm27, %zmm8
	vpermi2pd %zmm6, %zmm4, %zmm8
	vpermpd $0x44, %zmm6, %zmm2
	vmovapd %zmm27, %zmm9
	vpermi2pd %zmm4, %zmm6, %zmm9
	vfmadd231pd %zmm8, %zmm25, %zmm0
	vfmadd231pd %zmm9, %zmm26, %zmm2
	# Pair 2: (zmm5=re_hi, zmm7=im_hi) → (zmm4=new_re_hi, zmm6=new_im_hi)
	vpermpd $0x44, %zmm5, %zmm4
	vmovapd %zmm27, %zmm8
	vpermi2pd %zmm7, %zmm5, %zmm8
	vpermpd $0x44, %zmm7, %zmm6
	vmovapd %zmm27, %zmm9
	vpermi2pd %zmm5, %zmm7, %zmm9
	vfmadd231pd %zmm8, %zmm25, %zmm4
	vfmadd231pd %zmm9, %zmm26, %zmm6
	# Data: re_lo=zmm0, re_hi=zmm4, im_lo=zmm2, im_hi=zmm6
	# ── Size-8 (halfnn=4) butterfly across the two blocks ──
	vshuff64x2 $0x44, %zmm4, %zmm0, %zmm1
	vshuff64x2 $0xEE, %zmm4, %zmm0, %zmm3
	vshuff64x2 $0x44, %zmm6, %zmm2, %zmm5
	vshuff64x2 $0xEE, %zmm6, %zmm2, %zmm7
	vmulpd       %zmm3, %zmm28, %zmm8
	vfnmadd231pd %zmm7, %zmm29, %zmm8
	vmulpd       %zmm3, %zmm29, %zmm9
	vfmadd231pd  %zmm7, %zmm28, %zmm9
	vsubpd %zmm8, %zmm1, %zmm3
	vsubpd %zmm9, %zmm5, %zmm7
	vaddpd %zmm8, %zmm1, %zmm1
	vaddpd %zmm9, %zmm5, %zmm5
	vshuff64x2 $0x44, %zmm3, %zmm1, %zmm0
	vshuff64x2 $0xEE, %zmm3, %zmm1, %zmm4
	vshuff64x2 $0x44, %zmm7, %zmm5, %zmm2
	vshuff64x2 $0xEE, %zmm7, %zmm5, %zmm6
	# Store
	vmovapd %zmm0,    (%rdi,%rax,8)
	vmovapd %zmm4, 64(%rdi,%rax,8)
	vmovapd %zmm2,    (%rsi,%rax,8)
	vmovapd %zmm6, 64(%rsi,%rax,8)
	addq $16, %rax
	cmpq %r9, %rax
	jb fft_fused_248_loop
# ── Option C: fused {halfnn=8, halfnn=16, halfnn=32} pass ───────────────────
# Trig table layout relative to %r8 (base of runtime trig table):
#   [+0  ..+63 ] halfnn=4 entry (unused; fftsize8loop uses hardcoded consts)
#   [+64 ..+191] halfnn=8  entry: W8_cos at +64, W8_sin at +128
#   [+192..+447] halfnn=16 entry: W16[0]_cos at +192, W16[0]_sin at +256
#                                 W16[8]_cos at +320, W16[8]_sin at +384
#   [+448..+959] halfnn=32 entry: W32[0]_cos/sin at +448/+512
#                                 W32[8]_cos/sin at +576/+640
#                                 W32[16]_cos/sin at +704/+768
#                                 W32[24]_cos/sin at +832/+896
#   [+960..    ] halfnn=64+ entries → start of general loop
#
# Super-block of 64 doubles [A=+0, B=+8, C=+16, D=+24, E=+32, F=+40, G=+48, H=+56]
#   (each letter = 8 doubles = 1 ZMM)
#   Step 1 (halfnn=8):  A<->B, C<->D, E<->F, G<->H  (all use same W8)
#   Step 2 (halfnn=16): A'<->C', B'<->D', E'<->G', F'<->H'
#                       A'/C'/E'/G' use W16[0]; B'/D'/F'/H' use W16[8]
#   Step 3 (halfnn=32): A''<->E'', B''<->F'', C''<->G'', D''<->H''
#                       A''/E'' W32[0]; B''/F'' W32[8]; C''/G'' W32[16]; D''/H'' W32[24]
#
# Register map (all 32 ZMMs):
#   zmm0-7   = re data A..H  |  zmm8-15  = im data A..H
#   zmm16-17 = W8 cos/sin
#   zmm18-21 = W16[0]/W16[8] cos/sin
#   zmm22-29 = W32[0]/[8]/[16]/[24] cos/sin
#   zmm30-31 = temp re2/im2
	cmpq	$64,%r9
	jb	fft_fallback_to_8
	vmovapd  64(%r8),%zmm16		/* W8_cos */
	vmovapd 128(%r8),%zmm17		/* W8_sin */
	vmovapd 192(%r8),%zmm18		/* W16[0]_cos */
	vmovapd 256(%r8),%zmm19		/* W16[0]_sin */
	vmovapd 320(%r8),%zmm20		/* W16[8]_cos */
	vmovapd 384(%r8),%zmm21		/* W16[8]_sin */
	vmovapd 448(%r8),%zmm22		/* W32[0]_cos */
	vmovapd 512(%r8),%zmm23		/* W32[0]_sin */
	vmovapd 576(%r8),%zmm24		/* W32[8]_cos */
	vmovapd 640(%r8),%zmm25		/* W32[8]_sin */
	vmovapd 704(%r8),%zmm26		/* W32[16]_cos */
	vmovapd 768(%r8),%zmm27		/* W32[16]_sin */
	vmovapd 832(%r8),%zmm28		/* W32[24]_cos */
	vmovapd 896(%r8),%zmm29		/* W32[24]_sin */
	movq	$0,%rbx			/* rbx: super-block base (doubles) */
.p2align 4
fft_r4_8_16_32_loop:
	vmovapd    (%rdi,%rbx,8),%zmm0		/* A_re */
	vmovapd  64(%rdi,%rbx,8),%zmm1		/* B_re */
	vmovapd 128(%rdi,%rbx,8),%zmm2		/* C_re */
	vmovapd 192(%rdi,%rbx,8),%zmm3		/* D_re */
	vmovapd 256(%rdi,%rbx,8),%zmm4		/* E_re */
	vmovapd 320(%rdi,%rbx,8),%zmm5		/* F_re */
	vmovapd 384(%rdi,%rbx,8),%zmm6		/* G_re */
	vmovapd 448(%rdi,%rbx,8),%zmm7		/* H_re */
	vmovapd    (%rsi,%rbx,8),%zmm8		/* A_im */
	vmovapd  64(%rsi,%rbx,8),%zmm9		/* B_im */
	vmovapd 128(%rsi,%rbx,8),%zmm10	/* C_im */
	vmovapd 192(%rsi,%rbx,8),%zmm11	/* D_im */
	vmovapd 256(%rsi,%rbx,8),%zmm12	/* E_im */
	vmovapd 320(%rsi,%rbx,8),%zmm13	/* F_im */
	vmovapd 384(%rsi,%rbx,8),%zmm14	/* G_im */
	vmovapd 448(%rsi,%rbx,8),%zmm15	/* H_im */
	# halfnn=8: A vs B
	vmulpd	%zmm1,%zmm16,%zmm30
	vmulpd	%zmm1,%zmm17,%zmm31
	vfnmadd231pd %zmm9,%zmm17,%zmm30
	vfmadd231pd  %zmm9,%zmm16,%zmm31
	vsubpd	%zmm30,%zmm0,%zmm1
	vsubpd	%zmm31,%zmm8,%zmm9
	vaddpd	%zmm30,%zmm0,%zmm0
	vaddpd	%zmm31,%zmm8,%zmm8
	# halfnn=8: C vs D
	vmulpd	%zmm3,%zmm16,%zmm30
	vmulpd	%zmm3,%zmm17,%zmm31
	vfnmadd231pd %zmm11,%zmm17,%zmm30
	vfmadd231pd  %zmm11,%zmm16,%zmm31
	vsubpd	%zmm30,%zmm2,%zmm3
	vsubpd	%zmm31,%zmm10,%zmm11
	vaddpd	%zmm30,%zmm2,%zmm2
	vaddpd	%zmm31,%zmm10,%zmm10
	# halfnn=8: E vs F
	vmulpd	%zmm5,%zmm16,%zmm30
	vmulpd	%zmm5,%zmm17,%zmm31
	vfnmadd231pd %zmm13,%zmm17,%zmm30
	vfmadd231pd  %zmm13,%zmm16,%zmm31
	vsubpd	%zmm30,%zmm4,%zmm5
	vsubpd	%zmm31,%zmm12,%zmm13
	vaddpd	%zmm30,%zmm4,%zmm4
	vaddpd	%zmm31,%zmm12,%zmm12
	# halfnn=8: G vs H
	vmulpd	%zmm7,%zmm16,%zmm30
	vmulpd	%zmm7,%zmm17,%zmm31
	vfnmadd231pd %zmm15,%zmm17,%zmm30
	vfmadd231pd  %zmm15,%zmm16,%zmm31
	vsubpd	%zmm30,%zmm6,%zmm7
	vsubpd	%zmm31,%zmm14,%zmm15
	vaddpd	%zmm30,%zmm6,%zmm6
	vaddpd	%zmm31,%zmm14,%zmm14
	# halfnn=16: A' vs C'  (W16[0])
	vmulpd	%zmm2,%zmm18,%zmm30
	vmulpd	%zmm2,%zmm19,%zmm31
	vfnmadd231pd %zmm10,%zmm19,%zmm30
	vfmadd231pd  %zmm10,%zmm18,%zmm31
	vsubpd	%zmm30,%zmm0,%zmm2
	vsubpd	%zmm31,%zmm8,%zmm10
	vaddpd	%zmm30,%zmm0,%zmm0
	vaddpd	%zmm31,%zmm8,%zmm8
	# halfnn=16: B' vs D'  (W16[8])
	vmulpd	%zmm3,%zmm20,%zmm30
	vmulpd	%zmm3,%zmm21,%zmm31
	vfnmadd231pd %zmm11,%zmm21,%zmm30
	vfmadd231pd  %zmm11,%zmm20,%zmm31
	vsubpd	%zmm30,%zmm1,%zmm3
	vsubpd	%zmm31,%zmm9,%zmm11
	vaddpd	%zmm30,%zmm1,%zmm1
	vaddpd	%zmm31,%zmm9,%zmm9
	# halfnn=16: E' vs G'  (W16[0])
	vmulpd	%zmm6,%zmm18,%zmm30
	vmulpd	%zmm6,%zmm19,%zmm31
	vfnmadd231pd %zmm14,%zmm19,%zmm30
	vfmadd231pd  %zmm14,%zmm18,%zmm31
	vsubpd	%zmm30,%zmm4,%zmm6
	vsubpd	%zmm31,%zmm12,%zmm14
	vaddpd	%zmm30,%zmm4,%zmm4
	vaddpd	%zmm31,%zmm12,%zmm12
	# halfnn=16: F' vs H'  (W16[8])
	vmulpd	%zmm7,%zmm20,%zmm30
	vmulpd	%zmm7,%zmm21,%zmm31
	vfnmadd231pd %zmm15,%zmm21,%zmm30
	vfmadd231pd  %zmm15,%zmm20,%zmm31
	vsubpd	%zmm30,%zmm5,%zmm7
	vsubpd	%zmm31,%zmm13,%zmm15
	vaddpd	%zmm30,%zmm5,%zmm5
	vaddpd	%zmm31,%zmm13,%zmm13
	# halfnn=32: A'' vs E''  (W32[0])
	vmulpd	%zmm4,%zmm22,%zmm30
	vmulpd	%zmm4,%zmm23,%zmm31
	vfnmadd231pd %zmm12,%zmm23,%zmm30
	vfmadd231pd  %zmm12,%zmm22,%zmm31
	vsubpd	%zmm30,%zmm0,%zmm4
	vsubpd	%zmm31,%zmm8,%zmm12
	vaddpd	%zmm30,%zmm0,%zmm0
	vaddpd	%zmm31,%zmm8,%zmm8
	# halfnn=32: B'' vs F''  (W32[8])
	vmulpd	%zmm5,%zmm24,%zmm30
	vmulpd	%zmm5,%zmm25,%zmm31
	vfnmadd231pd %zmm13,%zmm25,%zmm30
	vfmadd231pd  %zmm13,%zmm24,%zmm31
	vsubpd	%zmm30,%zmm1,%zmm5
	vsubpd	%zmm31,%zmm9,%zmm13
	vaddpd	%zmm30,%zmm1,%zmm1
	vaddpd	%zmm31,%zmm9,%zmm9
	# halfnn=32: C'' vs G''  (W32[16])
	vmulpd	%zmm6,%zmm26,%zmm30
	vmulpd	%zmm6,%zmm27,%zmm31
	vfnmadd231pd %zmm14,%zmm27,%zmm30
	vfmadd231pd  %zmm14,%zmm26,%zmm31
	vsubpd	%zmm30,%zmm2,%zmm6
	vsubpd	%zmm31,%zmm10,%zmm14
	vaddpd	%zmm30,%zmm2,%zmm2
	vaddpd	%zmm31,%zmm10,%zmm10
	# halfnn=32: D'' vs H''  (W32[24])
	vmulpd	%zmm7,%zmm28,%zmm30
	vmulpd	%zmm7,%zmm29,%zmm31
	vfnmadd231pd %zmm15,%zmm29,%zmm30
	vfmadd231pd  %zmm15,%zmm28,%zmm31
	vsubpd	%zmm30,%zmm3,%zmm7
	vsubpd	%zmm31,%zmm11,%zmm15
	vaddpd	%zmm30,%zmm3,%zmm3
	vaddpd	%zmm31,%zmm11,%zmm11
	vmovapd %zmm0,   (%rdi,%rbx,8)
	vmovapd %zmm1, 64(%rdi,%rbx,8)
	vmovapd %zmm2,128(%rdi,%rbx,8)
	vmovapd %zmm3,192(%rdi,%rbx,8)
	vmovapd %zmm4,256(%rdi,%rbx,8)
	vmovapd %zmm5,320(%rdi,%rbx,8)
	vmovapd %zmm6,384(%rdi,%rbx,8)
	vmovapd %zmm7,448(%rdi,%rbx,8)
	vmovapd %zmm8,   (%rsi,%rbx,8)
	vmovapd %zmm9, 64(%rsi,%rbx,8)
	vmovapd %zmm10,128(%rsi,%rbx,8)
	vmovapd %zmm11,192(%rsi,%rbx,8)
	vmovapd %zmm12,256(%rsi,%rbx,8)
	vmovapd %zmm13,320(%rsi,%rbx,8)
	vmovapd %zmm14,384(%rsi,%rbx,8)
	vmovapd %zmm15,448(%rsi,%rbx,8)
	addq	$64,%rbx
	cmpq	%r9,%rbx
	jb	fft_r4_8_16_32_loop
	# Skip halfnn=4,8,16,32 trig entries (64+128+256+512=960 bytes)
	leaq	960(%r8),%rdx
	movq	$64,%rax
	# Option D: fuse halfnn=64+128 into a single pass over 4-column groups.
	# Each super-block covers 256 elements (4 columns × 64 elements).
	# Trig values repeat across super-blocks, so we reset trig pointers
	# each super-block and advance the data offset by 2048 bytes.
	cmpq	$256,%r9
	jb	ffthalfnnloop
	je	fft_d64_128_twist_setup	/* ns4==256: fuse with final twist */
	# ns4 > 256: Option D with outer super-block loop
	leaq	1024(%rdx),%r10		/* r10 = W128 lo base */
	leaq	1024(%r10),%r11		/* r11 = W128 hi base */
	movq	%rdx,%r12		/* r12 = saved W64 trig base */
	movq	%r10,%r13		/* r13 = saved W128 lo base */
	movq	%r11,%r14		/* r14 = saved W128 hi base */
	movq	$0,%rbx			/* rbx = data byte offset (outer+inner) */
	leaq	(,%r9,8),%rax		/* rax = ns4 * 8 = total data bytes */
.p2align 4
fft_d64_128_outer:
	movq	%r12,%rdx		/* reset W64 trig */
	movq	%r13,%r10		/* reset W128 lo trig */
	movq	%r14,%r11		/* reset W128 hi trig */
	leaq	512(%rbx),%rcx		/* rcx = end of this super-block inner */
.p2align 4
fft_d64_128_loop:
	# Load 4 columns of re data
	vmovapd	(%rdi,%rbx),%zmm0		/* re[i] */
	vmovapd	512(%rdi,%rbx),%zmm1		/* re[i+64] */
	vmovapd	1024(%rdi,%rbx),%zmm2		/* re[i+128] */
	vmovapd	1536(%rdi,%rbx),%zmm3		/* re[i+192] */
	# Load 4 columns of im data
	vmovapd	(%rsi,%rbx),%zmm4		/* im[i] */
	vmovapd	512(%rsi,%rbx),%zmm5		/* im[i+64] */
	vmovapd	1024(%rsi,%rbx),%zmm6		/* im[i+128] */
	vmovapd	1536(%rsi,%rbx),%zmm7		/* im[i+192] */
	# Load W64 trig (same for both halfnn=64 pairs)
	vmovapd	(%rdx),%zmm8			/* W64 cos */
	vmovapd	64(%rdx),%zmm9			/* W64 sin */
	# halfnn=64 butterfly A: zmm0<->zmm1 (re[i] <-> re[i+64])
	vmulpd	%zmm1,%zmm8,%zmm14
	vmulpd	%zmm1,%zmm9,%zmm15
	vfnmadd231pd %zmm5,%zmm9,%zmm14	/* re2 = re1*cos - im1*sin */
	vfmadd231pd  %zmm5,%zmm8,%zmm15	/* im2 = re1*sin + im1*cos */
	vsubpd	%zmm14,%zmm0,%zmm1		/* new re[i+64] = re[i] - re2 */
	vsubpd	%zmm15,%zmm4,%zmm5		/* new im[i+64] = im[i] - im2 */
	vaddpd	%zmm14,%zmm0,%zmm0		/* new re[i] = re[i] + re2 */
	vaddpd	%zmm15,%zmm4,%zmm4		/* new im[i] = im[i] + im2 */
	# halfnn=64 butterfly B: zmm2<->zmm3 (re[i+128] <-> re[i+192]), same W64
	vmulpd	%zmm3,%zmm8,%zmm14
	vmulpd	%zmm3,%zmm9,%zmm15
	vfnmadd231pd %zmm7,%zmm9,%zmm14
	vfmadd231pd  %zmm7,%zmm8,%zmm15
	vsubpd	%zmm14,%zmm2,%zmm3		/* new re[i+192] */
	vsubpd	%zmm15,%zmm6,%zmm7		/* new im[i+192] */
	vaddpd	%zmm14,%zmm2,%zmm2		/* new re[i+128] */
	vaddpd	%zmm15,%zmm6,%zmm6		/* new im[i+128] */
	# Load W128 trig (two different sets)
	vmovapd	(%r10),%zmm8			/* W128[i] cos */
	vmovapd	64(%r10),%zmm9			/* W128[i] sin */
	vmovapd	(%r11),%zmm10			/* W128[i+64] cos */
	vmovapd	64(%r11),%zmm11			/* W128[i+64] sin */
	# halfnn=128 butterfly A: zmm0<->zmm2 with W128[i]
	vmulpd	%zmm2,%zmm8,%zmm14
	vmulpd	%zmm2,%zmm9,%zmm15
	vfnmadd231pd %zmm6,%zmm9,%zmm14
	vfmadd231pd  %zmm6,%zmm8,%zmm15
	vsubpd	%zmm14,%zmm0,%zmm2		/* new re[i+128] */
	vsubpd	%zmm15,%zmm4,%zmm6		/* new im[i+128] */
	vaddpd	%zmm14,%zmm0,%zmm0		/* new re[i] */
	vaddpd	%zmm15,%zmm4,%zmm4		/* new im[i] */
	# halfnn=128 butterfly B: zmm1<->zmm3 with W128[i+64]
	vmulpd	%zmm3,%zmm10,%zmm14
	vmulpd	%zmm3,%zmm11,%zmm15
	vfnmadd231pd %zmm7,%zmm11,%zmm14
	vfmadd231pd  %zmm7,%zmm10,%zmm15
	vsubpd	%zmm14,%zmm1,%zmm3		/* new re[i+192] */
	vsubpd	%zmm15,%zmm5,%zmm7		/* new im[i+192] */
	vaddpd	%zmm14,%zmm1,%zmm1		/* new re[i+64] */
	vaddpd	%zmm15,%zmm5,%zmm5		/* new im[i+64] */
	# Store results
	vmovapd	%zmm0,(%rdi,%rbx)
	vmovapd	%zmm1,512(%rdi,%rbx)
	vmovapd	%zmm2,1024(%rdi,%rbx)
	vmovapd	%zmm3,1536(%rdi,%rbx)
	vmovapd	%zmm4,(%rsi,%rbx)
	vmovapd	%zmm5,512(%rsi,%rbx)
	vmovapd	%zmm6,1024(%rsi,%rbx)
	vmovapd	%zmm7,1536(%rsi,%rbx)
	# Advance inner
	addq	$64,%rbx		/* data: +8 doubles */
	leaq	128(%rdx),%rdx		/* W64 trig: +128 bytes */
	leaq	128(%r10),%r10		/* W128 lo: +128 bytes */
	leaq	128(%r11),%r11		/* W128 hi: +128 bytes */
	cmpq	%rcx,%rbx		/* end of this super-block? */
	jb	fft_d64_128_loop
	# Advance to next super-block: skip from sb+512 to sb+2048
	leaq	1536(%rbx),%rbx
	cmpq	%rax,%rbx		/* done with all super-blocks? */
	jb	fft_d64_128_outer
	# After fused pass: advance to halfnn=256 trig
	leaq	4032(%r8),%rdx		/* rdx = r8 + 960 + 1024 + 2048 */
	movq	$256,%rax		/* halfnn = 256 */
	cmpq	%r9,%rax
	jb	ffthalfnnloop		/* continue general loop for N >= 2048 */
	jmp	fftbeforefinal
fft_d64_128_twist_setup:
	# Fused Option D + final twist for N=1024 (ns4=256)
	# After butterflies, apply twist multiply in-register before storing.
	# Eliminates the separate fftfinalloop entirely.
	leaq	1024(%rdx),%r10		/* r10 = W128 base */
	leaq	1024(%r10),%r11		/* r11 = W128[i+64] base */
	leaq	4032(%r8),%rcx		/* rcx = twist trig base (r8+4032) */
	movq	$0,%rbx
.p2align 4
fft_d64_128_twist_loop:
	# Load 4 columns of data (identical to regular Option D)
	vmovapd	(%rdi,%rbx),%zmm0
	vmovapd	512(%rdi,%rbx),%zmm1
	vmovapd	1024(%rdi,%rbx),%zmm2
	vmovapd	1536(%rdi,%rbx),%zmm3
	vmovapd	(%rsi,%rbx),%zmm4
	vmovapd	512(%rsi,%rbx),%zmm5
	vmovapd	1024(%rsi,%rbx),%zmm6
	vmovapd	1536(%rsi,%rbx),%zmm7
	# W64 butterflies (identical to regular Option D)
	vmovapd	(%rdx),%zmm8
	vmovapd	64(%rdx),%zmm9
	vmulpd	%zmm1,%zmm8,%zmm14
	vmulpd	%zmm1,%zmm9,%zmm15
	vfnmadd231pd %zmm5,%zmm9,%zmm14
	vfmadd231pd  %zmm5,%zmm8,%zmm15
	vsubpd	%zmm14,%zmm0,%zmm1
	vsubpd	%zmm15,%zmm4,%zmm5
	vaddpd	%zmm14,%zmm0,%zmm0
	vaddpd	%zmm15,%zmm4,%zmm4
	vmulpd	%zmm3,%zmm8,%zmm14
	vmulpd	%zmm3,%zmm9,%zmm15
	vfnmadd231pd %zmm7,%zmm9,%zmm14
	vfmadd231pd  %zmm7,%zmm8,%zmm15
	vsubpd	%zmm14,%zmm2,%zmm3
	vsubpd	%zmm15,%zmm6,%zmm7
	vaddpd	%zmm14,%zmm2,%zmm2
	vaddpd	%zmm15,%zmm6,%zmm6
	# W128 butterflies (identical to regular Option D)
	vmovapd	(%r10),%zmm8
	vmovapd	64(%r10),%zmm9
	vmovapd	(%r11),%zmm10
	vmovapd	64(%r11),%zmm11
	vmulpd	%zmm2,%zmm8,%zmm14
	vmulpd	%zmm2,%zmm9,%zmm15
	vfnmadd231pd %zmm6,%zmm9,%zmm14
	vfmadd231pd  %zmm6,%zmm8,%zmm15
	vsubpd	%zmm14,%zmm0,%zmm2
	vsubpd	%zmm15,%zmm4,%zmm6
	vaddpd	%zmm14,%zmm0,%zmm0
	vaddpd	%zmm15,%zmm4,%zmm4
	vmulpd	%zmm3,%zmm10,%zmm14
	vmulpd	%zmm3,%zmm11,%zmm15
	vfnmadd231pd %zmm7,%zmm11,%zmm14
	vfmadd231pd  %zmm7,%zmm10,%zmm15
	vsubpd	%zmm14,%zmm1,%zmm3
	vsubpd	%zmm15,%zmm5,%zmm7
	vaddpd	%zmm14,%zmm1,%zmm1
	vaddpd	%zmm15,%zmm5,%zmm5
	# ── Apply final twist multiply to all 4 columns in-register ──
	# Column 0 (zmm0=re[i], zmm4=im[i]): twist at rcx
	vmovapd	(%rcx),%zmm8
	vmovapd	64(%rcx),%zmm9
	vmulpd	%zmm0,%zmm8,%zmm10
	vmulpd	%zmm0,%zmm9,%zmm11
	vfnmadd231pd %zmm4,%zmm9,%zmm10
	vfmadd231pd  %zmm4,%zmm8,%zmm11
	vmovapd	%zmm10,(%rdi,%rbx)
	vmovapd	%zmm11,(%rsi,%rbx)
	# Column 1 (zmm1=re[i+64], zmm5=im[i+64]): twist at rcx+1024
	vmovapd	1024(%rcx),%zmm8
	vmovapd	1088(%rcx),%zmm9
	vmulpd	%zmm1,%zmm8,%zmm10
	vmulpd	%zmm1,%zmm9,%zmm11
	vfnmadd231pd %zmm5,%zmm9,%zmm10
	vfmadd231pd  %zmm5,%zmm8,%zmm11
	vmovapd	%zmm10,512(%rdi,%rbx)
	vmovapd	%zmm11,512(%rsi,%rbx)
	# Column 2 (zmm2=re[i+128], zmm6=im[i+128]): twist at rcx+2048
	vmovapd	2048(%rcx),%zmm8
	vmovapd	2112(%rcx),%zmm9
	vmulpd	%zmm2,%zmm8,%zmm10
	vmulpd	%zmm2,%zmm9,%zmm11
	vfnmadd231pd %zmm6,%zmm9,%zmm10
	vfmadd231pd  %zmm6,%zmm8,%zmm11
	vmovapd	%zmm10,1024(%rdi,%rbx)
	vmovapd	%zmm11,1024(%rsi,%rbx)
	# Column 3 (zmm3=re[i+192], zmm7=im[i+192]): twist at rcx+3072
	vmovapd	3072(%rcx),%zmm8
	vmovapd	3136(%rcx),%zmm9
	vmulpd	%zmm3,%zmm8,%zmm10
	vmulpd	%zmm3,%zmm9,%zmm11
	vfnmadd231pd %zmm7,%zmm9,%zmm10
	vfmadd231pd  %zmm7,%zmm8,%zmm11
	vmovapd	%zmm10,1536(%rdi,%rbx)
	vmovapd	%zmm11,1536(%rsi,%rbx)
	# Advance
	addq	$64,%rbx
	leaq	128(%rdx),%rdx
	leaq	128(%r10),%r10
	leaq	128(%r11),%r11
	leaq	128(%rcx),%rcx
	cmpq	$512,%rbx
	jb	fft_d64_128_twist_loop
	jmp	fftend			/* skip fftfinalloop entirely */

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
.p2align 5
fftfinalloop:
	vmovapd	(%r10),%zmm0 /* re */
	vmovapd	(%r11),%zmm1 /* im */
	vmovapd (%rdx),%zmm2 /* cos */
	vmovapd 64(%rdx),%zmm3 /* sin */
	vmulpd %zmm0,%zmm2,%zmm4        /* re*cos */
	vmulpd %zmm0,%zmm3,%zmm5        /* re*sin */
	vfnmadd231pd %zmm1,%zmm3,%zmm4  /* re*cos - im*sin */
	vfmadd231pd %zmm1,%zmm2,%zmm5   /* re*sin + im*cos */
	vmovapd %zmm4,(%r10)
	vmovapd %zmm5,(%r11)
    	/* end of final loop */
    	leaq	64(%r10),%r10
    	leaq	64(%r11),%r11
    	leaq	128(%rdx),%rdx
	addq	$8,%rax
	cmpq	%r9,%rax
	jb fftfinalloop

	/* Restore registers */
fftend:
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
