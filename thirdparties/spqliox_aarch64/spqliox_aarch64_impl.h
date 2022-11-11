#pragma once
#include <xbyak_aarch64/xbyak_aarch64.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <memory>

namespace spqliox_aarch64 {

struct fft_code : Xbyak_aarch64::CodeGenerator {
    fft_code()
    {
        using namespace Xbyak_aarch64;
        // 退避
        // stp(x29, x30, AdrPreImm(sp, -16));
        // stp(x27, x28, AdrPreImm(sp, -16));
        // stp(x25, x26, AdrPreImm(sp, -16));
        // stp(x23, x24, AdrPreImm(sp, -16));
        // stp(x21, x22, AdrPreImm(sp, -16));
        stp(x19, x20, AdrPreImm(sp, -16));
        // sub(sp, sp, 64);
        // st1(v8.d2 - v11.d2, ptr(sp));
        // sub(sp, sp, 64);
        // st1(v12.d2 - v15.d2, ptr(sp));

        /*
    index  0    1    2    3    4      5      6
    reg    x0,  x1,  x2,  x3,  x4,    x5     x6
    var    dst, sit, send,bla, table, inout, negation_tbl
    */

        ld1r(v0.d2, ptr(x3));
        Label multiply;
        L(multiply);
        ld1(v1.d2, AdrPostImm(x1, 16));
        fmul(v1.d2, v0.d2, v1.d2);
        st1(v1.d2, AdrPostImm(x0, 16));
        cmp(x1, x2);
        bne(multiply);

        mov(x0, x4);
        mov(x1, x5);

        mov(x7, x0);
        mov(x0, x1);  // x0 = inout_direct

        // x7 is base of table
        ldr(x2, AdrPostImm(x7, 8));  // x2 is trig_tables.n
        // x4 is base of trig_tables
        ldr(x4, ptr(x7));

        // x2 is table N
        // x5 is N/4
        lsr(x5, x2, 2);
        lsl(x8, x2, 1);
        add(x1, x8, x0);  // x1 = 2N + inout_direct = &inout_direct[N/4]

        //    //size 2
        //    {
        //	//[1  1]
        //	//[1 -1]
        //	//     [1  1]
        //	//     [1 -1]
        //	for (int32_t block=0; block<ns4; block+=4) {
        //	    double* d0 = pre+block;
        //	    double* d1 = pim+block;
        //	    tmp0[0]=d0[0]; // d0[0] = v0.d2
        //	    tmp0[1]=d0[0];
        //	    tmp0[2]=d0[2]; // d0[2] = v2.d2
        //	    tmp0[3]=d0[2];
        //	    tmp1[0]=d0[1]; // d0[1] = v1.d2
        //	    tmp1[1]=-d0[1];
        //	    tmp1[2]=d0[3]; // d0[3] = v3.d2
        //	    tmp1[3]=-d0[3];
        //	    add4(d0,tmp0,tmp1); // d0[0:1] = v0.d2 + v1.d2 * v24.d2,
        // d[2:3] = v2.d2 + v3.d2 * v25.d2 	    tmp0[0]=d1[0]; // d1[0] =
        // v4.d2 	    tmp0[1]=d1[0]; 	    tmp0[2]=d1[2]; // d1[2] =
        // v6.d2 	    tmp0[3]=d1[2];
        //	    tmp1[0]=d1[1]; // d1[1] = v5.d2
        //	    tmp1[1]=-d1[1];
        //	    tmp1[2]=d1[3]; // d1[3] = v7.d2
        //	    tmp1[3]=-d1[3];
        //	    add4(d1,tmp0,tmp1); // d1[0:1] = v4.d2 + v5.d2 * v25.d2,
        // d[2:3] = v6.d2 + v7.d2 * v25.d2
        //	}
        //    }

        // x6 is base of negation table
        ld1(v24.d2, AdrPostImm(x6, 16));
        ld1(v25.d2, AdrPostImm(x6, 16));
        ld1(v26.d2, AdrPostImm(x6, 16));
        ld1(v27.d2, AdrPostImm(x6, 16));
        ld1(v28.d2, AdrPostImm(x6, 16));
        ld1(v29.d2, AdrPostImm(x6, 16));

        mov(x7, 0);
        mov(x9, x0);   // x9 = inout_direct
        mov(x10, x1);  // x10 = &inout_direct[N/4]
        Label fftsize2loop;
        L(fftsize2loop);
        ld1(v16.d2 - v17.d2, ptr(x9));
        ld1(v18.d2 - v19.d2, ptr(x10));
        dup(v0.d2, v16.d2[0]);
        dup(v2.d2, v16.d2[1]);
        dup(v1.d2, v17.d2[0]);
        dup(v3.d2, v17.d2[1]);
        dup(v4.d2, v18.d2[0]);
        dup(v6.d2, v18.d2[1]);
        dup(v5.d2, v19.d2[0]);
        dup(v7.d2, v19.d2[1]);
        fmla(v0.d2, v2.d2, v24.d2);
        fmla(v1.d2, v3.d2, v24.d2);
        st1(v0.d2 - v1.d2, AdrPostImm(x9, 32));
        fmla(v4.d2, v6.d2, v25.d2);
        fmla(v5.d2, v7.d2, v25.d2);
        st1(v4.d2 - v5.d2, AdrPostImm(x10, 32));
        add(x7, x7, 4);
        cmp(x7, x5);
        blt(fftsize2loop);

        //
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
        //  re[1:0] = v0.d2
        //  re[3:2] = v1.d2
        //  im[1:0] = v2.d2
        //  im[3:2] = v3.d2
        //  tmp1[1:0](v4.d2) = {v3.d2[1], v1.d2[0]}
        //  v0.d2(re[1:0]) = v0.d2 + tmp1[1:0]
        //  v5.d2 = tmp1[1:0] * v27.d2
        //  v1.d2(re[3:2]) = v0.d2 + v5.d2
        //  tmp3[1:0](v6.d2) = {v1.d2[1], v3.d2[0]}
        //  v7.d2 = tmp3[1:0] * v28.d2
        //  v8.d2 = tmp3[1:0] * v29.d2
        //  v2.d2(im[1:0]) = v2.d2 + v7.d2
        //  v3.d2(im[3:2]) = v2.d2 + v8.d2
        //
        //  re[3:0] = tmp0[3:0] + tmp1[3:0]
        //	for (int32_t block=0; block<ns4; block+=4) {
        //	    double* re = pre+block;
        //	    double* im = pim+block;
        //	    tmp0[0]=re[0];
        //	    tmp0[1]=re[1]; //re[0:1] = v0.d2
        //	    tmp0[2]=re[0];
        //	    tmp0[3]=re[1];
        //	    tmp1[0]=re[2];
        //	    tmp1[1]=im[3]; //v2.d2
        //	    tmp1[2]=-re[2];
        //	    tmp1[3]=-im[3]; //v3.d2
        //	    tmp2[0]=im[0];
        //	    tmp2[1]=im[1]; // v4.d2
        //	    tmp2[2]=im[0];
        //	    tmp2[3]=im[1]; // v5.d2
        //	    tmp3[0]=im[2];
        //	    tmp3[1]=-re[3]; // v6.d2
        //	    tmp3[2]=-im[2];
        //	    tmp3[3]=re[3]; // v7.d2
        //	    add4(re,tmp0,tmp1);
        //	    add4(im,tmp2,tmp3);
        //	}
        //    }
        mov(x7, 0);
        mov(x9, x0);   // x9 = inout_direct
        mov(x10, x1);  // x10 = &inout_direct[N/4]
        Label fftsize4loop;
        L(fftsize4loop);
        ld1(v0.d2 - v1.d2, ptr(x9));
        ld1(v2.d2 - v3.d2, ptr(x10));
        ins(v4.d2[1], v3.d2[1]);
        ins(v4.d2[0], v1.d2[0]);
        ins(v6.d2[1], v1.d2[1]);
        ins(v6.d2[0], v3.d2[0]);
        fmul(v5.d2, v4.d2, v27.d2);
        fadd(v1.d2, v0.d2, v5.d2);
        fadd(v0.d2, v0.d2, v4.d2);
        st1(v0.d2 - v1.d2, AdrPostImm(x9, 32));
        fmul(v7.d2, v6.d2, v28.d2);
        fmul(v16.d2, v6.d2, v29.d2);
        fadd(v3.d2, v2.d2, v16.d2);
        fadd(v2.d2, v2.d2, v7.d2);
        st1(v2.d2 - v3.d2, AdrPostImm(x10, 32));
        add(x7, x7, 4);
        cmp(x7, x5);
        blt(fftsize4loop);

        //
        ////general loop
        // const double* cur_tt = trig_tables; //x2
        // for (int32_t halfnn=4 //x7 ; halfnn<ns4; halfnn*=2) {
        //	int32_t nn = 2*halfnn; /x11
        //	for (int32_t block=0 //x12; block<ns4; block+=nn) {
        //	    for (int32_t off=0 //x13; off<halfnn; off+=4) {
        //		    double* re0 = pre + block + off; // x14 (re0
        //+= 0x20) 		    double* im0 = pim + block + off; // x15 (im0
        //+=
        // 0x20) 		    double* re1 = pre + block + halfnn + off; //
        // x16
        // (re1 += 0x20) 		    double* im1 = pim + block + halfnn +
        // off;
        // // x17 (im1
        //+= 0x20) 		    const double* tcs = cur_tt+2*off;  // x18
        //(tcs += 0x40) 		    const double* tsn = tcs+4; // (tcs +
        // 0x20) 		    dotp4(tmp0,re1,tcs); // re*cos
        //		    dotp4(tmp1,re1,tsn); // re*sin
        //		    dotp4(tmp2,im1,tcs); // im*cos
        //		    dotp4(tmp3,im1,tsn); // im*sin
        //		    sub4(tmp0,tmp0,tmp3); // re2
        //		    add4(tmp1,tmp1,tmp2); // im2
        //		    add4(tmp2,re0,tmp0); // re + re
        //		    add4(tmp3,im0,tmp1); // im + im
        //		    sub4(tmp0,re0,tmp0); // re - re
        //		    sub4(tmp1,im0,tmp1); // im - im
        //		    copy4(re0,tmp2);
        //		    copy4(im0,tmp3);
        //		    copy4(re1,tmp0);
        //		    copy4(im1,tmp1);
        //	    }
        //	}
        //	cur_tt += nn;
        //}

        mov(x2, x4);  // x4 is base of trig_tables
        mov(x7, 4);
        Label ffthalfnnloop;
        L(ffthalfnnloop);
        lsl(x11, x7, 1);
        Label fftblockloop;
        mov(x12, 0);
        lsl(x19, x7, 3);  // ptr(halfnn) = halfnn * 8
        L(fftblockloop);
        mov(x13, 0);
        lsl(x20, x12, 3);  // ptr(block) = block * 8
        add(x14, x0, x20);
        add(x15, x1, x20);
        add(x16, x14, x19);
        add(x17, x15, x19);
        mov(x18, x2);
        Label fftoffloop;
        L(fftoffloop);
        ld1(v0.d2 - v1.d2, ptr(x14));  // re0[1:0]
        ld1(v2.d2 - v3.d2, ptr(x15));  // im0[1:0]
        ld1(v4.d2 - v5.d2, ptr(x16));  // re1[1:0]
        ld1(v6.d2 - v7.d2, ptr(x17));  // im1[1:0]
        // avoid use of v8 - v15 (because of cost for stack ops)
        ld1(v20.d2 - v23.d2, AdrPostImm(x18, 64));  // tcs[1:0]
        fmul(v24.d2, v4.d2, v20.d2);                // tmp0[1:0]
        fmul(v25.d2, v5.d2, v21.d2);                // tmp0[3:2]
        fmul(v26.d2, v4.d2, v22.d2);                // tmp1[1:0]
        fmul(v27.d2, v5.d2, v23.d2);                // tmp1[3:2]
        fmul(v16.d2, v6.d2, v20.d2);                // tmp2[1:0]
        fmul(v17.d2, v7.d2, v21.d2);                // tmp2[3:2]
        fmul(v18.d2, v6.d2, v22.d2);                // tmp3[1:0]
        fmul(v19.d2, v7.d2, v23.d2);                // tmp3[3:2]
        fsub(v24.d2, v24.d2, v18.d2);
        fsub(v25.d2, v25.d2, v19.d2);  // re2
        fadd(v26.d2, v26.d2, v16.d2);
        fadd(v27.d2, v27.d2, v17.d2);  // im2
        fadd(v16.d2, v0.d2, v24.d2);
        fadd(v17.d2, v1.d2, v25.d2);  // re + re
        fadd(v18.d2, v2.d2, v26.d2);
        fadd(v19.d2, v3.d2, v27.d2);  // im + im
        fsub(v24.d2, v0.d2, v24.d2);
        fsub(v25.d2, v1.d2, v25.d2);  // re - re
        fsub(v26.d2, v2.d2, v26.d2);
        fsub(v27.d2, v3.d2, v27.d2);  // im - im
        st1(v16.d2 - v17.d2, AdrPostImm(x14, 32));
        st1(v18.d2 - v19.d2, AdrPostImm(x15, 32));
        st1(v24.d2 - v25.d2, AdrPostImm(x16, 32));
        st1(v26.d2 - v27.d2, AdrPostImm(x17, 32));
        add(x13, x13, 4);
        cmp(x13, x7);
        blt(fftoffloop);
        add(x12, x12, x11);
        cmp(x12, x5);
        blt(fftblockloop);
        lsl(x7, x7, 1);
        lsl(x20, x11, 3);
        add(x2, x2, x20);  // cur_tt += n
        cmp(x7, x5);
        blt(ffthalfnnloop);

        //    //multiply by omb^j
        // for (int32_t j=0 //x7; j<ns4; j+=4) {
        //	const double* r0 = cur_tt+2*j; //x2
        //	const double* r1 = r0+4;
        //	//(re*cos-im*sin) + i (im*cos+re*sin)
        //	double* d0 = pre+j; //x10 pre x0
        //	double* d1 = pim+j; //x11 pim x1
        //	dotp4(tmp0,d0,r0); //re*cos
        //	dotp4(tmp1,d1,r0); //im*cos
        //	dotp4(tmp2,d0,r1); //re*sin
        //	dotp4(tmp3,d1,r1); //im*sin
        //	sub4(d0,tmp0,tmp3);
        //	add4(d1,tmp1,tmp2);
        //    }
        //}

        mov(x7, 0);
        mov(x10, x0);
        mov(x11, x1);
        Label fftfinalloop;
        L(fftfinalloop);
        ld1(v0.d2 - v3.d2, AdrPostImm(x2, 64));  // r0[1:0]
        ld1(v4.d2 - v5.d2, ptr(x10));            // d0[1:0]
        ld1(v6.d2 - v7.d2, ptr(x11));            // d1[1:0]
        fmul(v16.d2, v4.d2, v0.d2);              // tmp0[1:0]
        fmul(v17.d2, v5.d2, v1.d2);              // tmp0[3:2]
        fmul(v18.d2, v6.d2, v0.d2);              // tmp1[1:0]
        fmul(v19.d2, v7.d2, v1.d2);              // tmp1[3:2]
        fmul(v20.d2, v4.d2, v2.d2);              // tmp2[1:0]
        fmul(v21.d2, v5.d2, v3.d2);              // tmp2[3:2]
        fmul(v22.d2, v6.d2, v2.d2);              // tmp3[1:0]
        fmul(v23.d2, v7.d2, v3.d2);              // tmp3[3:2]
        fsub(v4.d2, v16.d2, v22.d2);
        fsub(v5.d2, v17.d2, v23.d2);
        fadd(v6.d2, v18.d2, v20.d2);
        fadd(v7.d2, v19.d2, v21.d2);
        st1(v4.d2 - v5.d2, AdrPostImm(x10, 32));
        st1(v6.d2 - v7.d2, AdrPostImm(x11, 32));
        add(x7, x7, 4);
        cmp(x7, x5);
        blt(fftfinalloop);

        // 復帰
        // v8 - v15 must be preserved
        // ld1(v12.d2 - v15.d2, AdrPostImm(sp, 64));
        // ld1(v8.d2 - v11.d2, AdrPostImm(sp, 64));
        // x19 - x30 must be preserved
        ldp(x19, x20, AdrPostImm(sp, 16));
        // ldp(x21, x22, AdrPostImm(sp, 16));
        // ldp(x23, x24, AdrPostImm(sp, 16));
        // ldp(x25, x26, AdrPostImm(sp, 16));
        // ldp(x27, x28, AdrPostImm(sp, 16));
        // ldp(x29, x30, AdrPostImm(sp, 16));
        ret();
    }
};

struct ifft_code : Xbyak_aarch64::CodeGenerator {
    ifft_code()
    {
        /*
           index  0    1    2     3      4      5
           reg    x0,  x1,  x2,   x3,    x4,    x5
           var    dst  ait, aend, table, inout, neg_table
        */
        using namespace Xbyak_aarch64;
        // 退避
        // v8 - v15 must be preserved
        // x19 - x30 must be preserved
        // stp(x29, x30, AdrPreImm(sp, -16));
        // stp(x27, x28, AdrPreImm(sp, -16));
        // stp(x25, x26, AdrPreImm(sp, -16));
        // stp(x23, x24, AdrPreImm(sp, -16));
        // stp(x21, x22, AdrPreImm(sp, -16));
        stp(x19, x20, AdrPreImm(sp, -16));

        ldr(x0, AdrPostImm(x3, 8));  // x0 is trig_tables.n = N
        lsr(x1, x0, 2);              // x1 is N/4
        ldr(x2, ptr(x3));            // x2 is base of trig_tables
        lsl(x7, x0, 1);              // x7 is 2N
        add(x7, x7, x4);  // x7 = 2N + inout_direct = &inout_direct[N/4]

        //    //multiply by omega^j
        // for (int32_t j=0 //x8; j<ns4; j+=4) {
        //	const double* r0 = trig_tables+2*j; //v0.d2, v1.d2 (=
        // trig_tables + 0x40) 	const double* r1 = r0+4; //v2.d2, v3.d2
        //	//(re*cos-im*sin) + i (im*cos+re*sin)
        //	double* d0 = are+j; // v4.d2, v5.d2
        //	double* d1 = aim+j; // v6.d2, v7.d2
        //	dotp4(tmp0,d0,r0); //re*cos
        //	dotp4(tmp1,d1,r0); //im*cos
        //	dotp4(tmp2,d0,r1); //re*sin
        //	dotp4(tmp3,d1,r1); //im*sin
        //	sub4(d0,tmp0,tmp3);
        //	add4(d1,tmp1,tmp2);
        //}

        // x1 = N/4
        // x2 = trig_tables
        // x4 = inout_direct
        // x7 = &inout_direct[N/4]
        mov(x8, 0);
        mov(x9, x2);
        mov(x10, x4);
        mov(x11, x7);
        mov(x0, x4);
        Label fftstartloop;
        L(fftstartloop);
        ld1(v0.d2 - v3.d2, AdrPostImm(x9, 64));  // r0, r1
        ld1(v4.d2 - v5.d2, ptr(x10));            // d0
        ld1(v6.d2 - v7.d2, ptr(x11));            // d1
        // avoid use of v8 - v15 because of cost for stack ops
        fmul(v16.d2, v4.d2, v0.d2);  // tmp0[1:0]
        fmul(v17.d2, v5.d2, v1.d2);  // tmp0[3:2]
        fmul(v18.d2, v6.d2, v0.d2);  // tmp1[1:0]
        fmul(v19.d2, v7.d2, v1.d2);  // tmp1[3:2]
        fmul(v20.d2, v4.d2, v2.d2);  // tmp2[1:0]
        fmul(v21.d2, v5.d2, v3.d2);  // tmp2[3:2]
        fmul(v22.d2, v6.d2, v2.d2);  // tmp3[1:0]
        fmul(v23.d2, v7.d2, v3.d2);  // tmp3[3:2]
        fsub(v4.d2, v16.d2, v22.d2);
        fsub(v5.d2, v17.d2, v23.d2);
        fadd(v6.d2, v18.d2, v20.d2);
        fadd(v7.d2, v19.d2, v21.d2);
        st1(v4.d2 - v5.d2, AdrPostImm(x10, 32));
        st1(v6.d2 - v7.d2, AdrPostImm(x11, 32));
        add(x8, x8, 4);
        cmp(x8, x1);
        blt(fftstartloop);

        /*
            const double* cur_tt = trig_tables; // x8
            for (int32_t nn=ns4 //x9; nn>=8; nn/=2) {
                int32_t halfnn = nn/2; //x10
                cur_tt += 2*nn;
                for (int32_t block=0 //x11; block<ns4; block+=nn) {
                    for (int32_t off=0 //x12; off<halfnn; off+=4) {
                        double* d00 = are + block + off; // x13
                        double* d01 = aim + block + off; // x14
                        double* d10 = are + block + halfnn + off; // x15
                        double* d11 = aim + block + halfnn + off; // x16
                        add4(tmp0,d00,d10); // re + re
                        add4(tmp1,d01,d11); // im + im
                        sub4(tmp2,d00,d10); // re - re
                        sub4(tmp3,d01,d11); // im - im
                        copy4(d00,tmp0);
                        copy4(d01,tmp1);
                        const double* r0 = cur_tt+2*off; // x17
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

        // x1 = N/4
        // x2 = trig_tables
        // x4 = inout_direct
        // x7 = &inout_direct[N/4]
        mov(x8, x2);
        mov(x9, x1);
        Label fftnnloop, fftblockloop, fftoffloop;
        L(fftnnloop);
        lsr(x10, x9, 1);   // x10 = nn/2
        lsl(x11, x9, 4);   // x11 = 2*nn *8
        add(x8, x8, x11);  // cur_tt += 2*nn*8
        mov(x11, 0);
        lsl(x19, x10, 3);  // x19 = ptr(hlafnn) = halfnn * 8
        L(fftblockloop);
        mov(x12, 0);
        lsl(x20, x11, 3);  // x20 = ptr(block) = block * 8
        add(x13, x4, x20);
        add(x14, x7, x20);
        add(x15, x13, x19);
        add(x16, x14, x19);
        mov(x17, x8);
        L(fftoffloop);
        ld1(v0.d2 - v1.d2, ptr(x13));  // d00
        ld1(v2.d2 - v3.d2, ptr(x14));  // d01
        ld1(v4.d2 - v5.d2, ptr(x15));  // d10
        ld1(v6.d2 - v7.d2, ptr(x16));  // d11
        // avoid use of v8 - v15 (because of cost for stack ops)
        ld1(v16.d2 - v19.d2, AdrPostImm(x17, 64));  // r0, r1
        fadd(v20.d2, v0.d2, v4.d2);                 // tmp0[1:0]
        fadd(v21.d2, v1.d2, v5.d2);                 // tmp0[3:2]
        fadd(v22.d2, v2.d2, v6.d2);                 // tmp1[1:0]
        fadd(v23.d2, v3.d2, v7.d2);                 // tmp1[3:2]
        fsub(v24.d2, v0.d2, v4.d2);                 // tmp2[1:0]
        fsub(v25.d2, v1.d2, v5.d2);                 // tmp2[3:2]
        fsub(v26.d2, v2.d2, v6.d2);                 // tmp3[1:0]
        fsub(v27.d2, v3.d2, v7.d2);                 // tmp3[3:2]
        st1(v20.d2 - v21.d2, AdrPostImm(x13, 32));
        st1(v22.d2 - v23.d2, AdrPostImm(x14, 32));
        fmul(v20.d2, v24.d2, v16.d2);
        fmul(v21.d2, v25.d2, v17.d2);
        fmul(v22.d2, v26.d2, v18.d2);
        fmul(v23.d2, v27.d2, v19.d2);
        fsub(v4.d2, v20.d2, v22.d2);
        fsub(v5.d2, v21.d2, v23.d2);
        st1(v4.d2 - v5.d2, AdrPostImm(x15, 32));
        fmul(v20.d2, v24.d2, v18.d2);
        fmul(v21.d2, v25.d2, v19.d2);
        fmul(v22.d2, v26.d2, v16.d2);
        fmul(v23.d2, v27.d2, v17.d2);
        fadd(v6.d2, v20.d2, v22.d2);
        fadd(v7.d2, v21.d2, v23.d2);
        st1(v6.d2 - v7.d2, AdrPostImm(x16, 32));
        add(x12, x12, 4);
        cmp(x12, x10);
        blt(fftoffloop);
        add(x11, x11, x9);
        cmp(x11, x1);
        blt(fftblockloop);
        mov(x9, x10);
        cmp(x9, 8);
        bge(fftnnloop);

        /*
            //size 4 loop
            {
                for (int32_t block=0 //x8; block<ns4; block+=4) {
                    double* d0 = are+block; v0.d2, v1.d2
                    double* d1 = aim+block; v2.d2, v3.d2
                    tmp0[0]=d0[0];
                    tmp0[1]=d0[1]; // v4.d2 = v0.d2
                    tmp0[2]=d0[0];
                    tmp0[3]=-d1[1]; // v5.d2 = {v2.d2[1], v0.d2[0]} * v25.d2
                    tmp1[0]=d0[2];
                    tmp1[1]=d0[3]; // v6.d2 = v1.d2
                    tmp1[2]=-d0[2];
                    tmp1[3]=d1[3]; // v7.d2 = {v3.d2[1], v1.d2[0]} * v27.d2
                    tmp2[0]=d1[0];
                    tmp2[1]=d1[1]; // v20.d2 = v2.d2
                    tmp2[2]=d1[0];
                    tmp2[3]=d0[1]; // v21.d2 = {v0.d2[1], v2.d2[0]}
                    tmp3[0]=d1[2];
                    tmp3[1]=d1[3]; // v22.d2 = v3.d2
                    tmp3[2]=-d1[2];
                    tmp3[3]=-d0[3]; // v23.d2 = {v1.d2[1], v3.d2[0]}
                    add4(d0,tmp0,tmp1); // v16.d2 = v0.d2 + v1.d2, v17.d2 =
           v5.d2 + v7.d2 add4(d1,tmp2,tmp3); // v18.d2 = v2.d2 + v3.d2, v19.d2 =
           v20.d2 - v11.d2
                }
            }
                    tmp0[0] = re[0] + re[2];
                    tmp0[1] = re[1] + re[3];
                    tmp0[2] = re[0] - re[2];
                    tmp0[3] = im[3] - im[1];

                    tmp1[0] = im[0] + im[2];
                    tmp1[1] = im[1] + im[3];
                    tmp1[2] = im[0] + -im[2];
                    tmp1[3] = re[1] + -re[3];

        */
        // x5 is base of negation table
        ld1(v24.d2, AdrPostImm(x5, 16));
        ld1(v25.d2, AdrPostImm(x5, 16));
        ld1(v26.d2, AdrPostImm(x5, 16));
        ld1(v27.d2, AdrPostImm(x5, 16));
        ld1(v28.d2, AdrPostImm(x5, 16));
        ld1(v29.d2, AdrPostImm(x5, 16));

        // x1 = N/4
        // x2 = trig_tables
        // x4 = inout_direct
        // x7 = &inout_direct[N/4]
        mov(x8, 0);
        mov(x9, x4);
        mov(x10, x7);
        Label ifftsize4loop;
        L(ifftsize4loop);
        ld1(v0.d2 - v1.d2, ptr(x9));
        ld1(v2.d2 - v3.d2, ptr(x10));
        ins(v5.d2[0], v0.d2[0]);
        ins(v5.d2[1], v2.d2[1]);
        ins(v7.d2[0], v1.d2[0]);
        ins(v7.d2[1], v3.d2[1]);
        ins(v20.d2[0], v2.d2[0]);
        ins(v20.d2[1], v0.d2[1]);
        ins(v11.d2[0], v3.d2[0]);
        ins(v11.d2[1], v1.d2[1]);
        fmul(v5.d2, v5.d2, v25.d2);
        fmul(v7.d2, v7.d2, v27.d2);
        fadd(v16.d2, v0.d2, v1.d2);
        fadd(v17.d2, v5.d2, v7.d2);
        fadd(v18.d2, v2.d2, v3.d2);
        fsub(v19.d2, v20.d2, v11.d2);
        st1(v16.d2 - v17.d2, AdrPostImm(x9, 32));
        st1(v18.d2 - v19.d2, AdrPostImm(x10, 32));
        add(x8, x8, 4);
        cmp(x8, x1);
        blt(ifftsize4loop);

        /*
            //size 2
                for (int32_t block=0; block<ns4; block+=4) {
                    double* d0 = are+block; v0.d2, v1.d2
                    double* d1 = aim+block; v2.d2, v3.d2
                    tmp0[0]=d0[0];
                    tmp0[1]=d0[0]; // v16.d2 = v0.d2[0]
                    tmp0[2]=d0[2];
                    tmp0[3]=d0[2]; // v17.d2 = v1.d2[0]
                    tmp1[0]=d0[1];
                    tmp1[1]=-d0[1]; // v18.d2 = v0.d2[1] * v25.d2
                    tmp1[2]=d0[3];
                    tmp1[3]=-d0[3]; // v19.d2 = v1.d2[1] * v25.d2
                    add4(d0,tmp0,tmp1); v20.d2 = v16.d2 + v18.d2, v21.d2 =
           v17.d2 + v19.d2 tmp0[0]=d1[0]; tmp0[1]=d1[0]; // v22.d2 = v2.d2[0]
                    tmp0[2]=d1[2];
                    tmp0[3]=d1[2]; // v23.d2 = v3.d2[0]
                    tmp1[0]=d1[1];
                    tmp1[1]=-d1[1]; // v4.d2 = v2.d2[1] * v25.d2
                    tmp1[2]=d1[3];
                    tmp1[3]=-d1[3]; // v5.d2 = v3.d2[1] * v25.d2
                    add4(d1,tmp0,tmp1); // v6.d2 = v22.d2 + v4.d2, v7.d2 =
           v23.d2 + v5.d2
                }

            tmp0[0] = re[0] + re[1];
            tmp0[1] = re[0] - re[1];
            tmp0[2] = re[2] + re[3];
            tmp0[3] = re[2] - re[3];
            copy4(re, tmp0);

            tmp1[0] = im[0] + im[1];
            tmp1[1] = im[0] - im[1];
            tmp1[2] = im[2] + im[3];
            tmp1[3] = im[2] - im[3];
            copy4(im, tmp1);
        */
        // x1 = N/4
        // x2 = trig_tables
        // x4 = inout_direct
        // x7 = &inout_direct[N/4]
        mov(x8, 0);
        mov(x9, x4);
        mov(x10, x7);
        Label ifftsize2loop;
        L(ifftsize2loop);
        ld1(v0.d2 - v1.d2, ptr(x9));
        ld1(v2.d2 - v3.d2, ptr(x10));
        dup(v16.d2, v0.d2[0]);
        dup(v17.d2, v1.d2[0]);
        dup(v18.d2, v0.d2[1]);
        dup(v19.d2, v1.d2[1]);
        dup(v22.d2, v2.d2[0]);
        dup(v23.d2, v3.d2[0]);
        dup(v4.d2, v2.d2[1]);
        dup(v5.d2, v3.d2[1]);
        fmul(v18.d2, v18.d2, v25.d2);
        fmul(v19.d2, v19.d2, v25.d2);
        fmul(v4.d2, v4.d2, v25.d2);
        fmul(v5.d2, v5.d2, v25.d2);
        fadd(v20.d2, v16.d2, v18.d2);
        fadd(v21.d2, v17.d2, v19.d2);
        st1(v20.d2 - v21.d2, AdrPostImm(x9, 32));
        fadd(v6.d2, v22.d2, v4.d2);
        fadd(v7.d2, v23.d2, v5.d2);
        st1(v6.d2 - v7.d2, AdrPostImm(x10, 32));
        add(x8, x8, 4);
        cmp(x8, x1);
        blt(ifftsize2loop);

        //// 復帰
        // v8 - v15 must be preserved
        // x19 - x30 must be preserved
        ldp(x19, x20, AdrPostImm(sp, 16));
        // ldp(x21, x22, AdrPostImm(sp, 16));
        // ldp(x23, x24, AdrPostImm(sp, 16));
        // ldp(x25, x26, AdrPostImm(sp, 16));
        // ldp(x27, x28, AdrPostImm(sp, 16));
        // ldp(x29, x30, AdrPostImm(sp, 16));
        ret();
    }
};
}  // namespace spqliox_aarch64